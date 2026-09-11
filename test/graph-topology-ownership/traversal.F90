! One source builds against the original pointer interface (-DBASELINE) or
! the read-only fibre interface. Construction is outside every measurement.
module traversal_sums
  use iso_c_binding, only : c_int64_t
  implicit none
  integer(c_int64_t) :: checksum
  type :: fibre_sum
     integer(c_int64_t) :: value = 0
  end type fibre_sum
contains
  subroutine sum_members(members)
    integer, intent(in) :: members(:)
    integer :: index
    do index = 1, size(members)
       checksum = checksum + members(index)
    end do
  end subroutine sum_members

  subroutine sum_context_members(members, context)
    integer, intent(in) :: members(:)
    class(*), intent(inout) :: context
    integer :: index
    select type(context)
    type is(fibre_sum)
       do index = 1, size(members)
          context % value = context % value + members(index)
       end do
    class default
       error stop 'the fibre sum context is required'
    end select
  end subroutine sum_context_members

  subroutine sum_incoming(offsets, indices, sources)
    integer, intent(in) :: offsets(:), indices(:), sources(:)
    integer :: vertex, position, edge_index
    do vertex = 1, size(offsets)-1
       do position = offsets(vertex), offsets(vertex+1)-1
          edge_index = indices(position)
          checksum = checksum + edge_index + sources(edge_index) + vertex
       end do
    end do
  end subroutine sum_incoming
end module traversal_sums

program topology_traversal

  use iso_c_binding, only : c_int64_t
  use iso_fortran_env, only : real64
  use graph_fractal, only : graph
  use map_set, only : set_map
  use map_set_representation, only : counted_set_representation
  use relation_binary, only : csr_relation
#ifndef BASELINE
  use relation_binary, only : integer_fibre
#endif
  use view_directed_stored, only : stored_directed_graph
  use traversal_sums, only : checksum, sum_members, sum_incoming, fibre_sum, sum_context_members

  implicit none

  interface
     subroutine allocation_begin() bind(C)
     end subroutine allocation_begin
     function allocation_end() result(num_allocations) bind(C)
       import c_int64_t
       integer(c_int64_t) :: num_allocations
     end function allocation_end
  end interface

  integer, parameter :: degree = 6
  integer :: num_vertices, num_repetitions, i, j, e, repetition
  integer, allocatable :: table(:,:), tails(:), heads(:), owned(:)
  integer(c_int64_t) :: expected_sum, num_allocations, graph_sum
  integer(c_int64_t) :: initial_count, final_count, count_rate
  real(real64) :: seconds
  character(len=32) :: argument
  type(graph) :: domain
  type(set_map) :: sets
  type(csr_relation), target :: relation
  type(stored_directed_graph) :: incidence
  type(fibre_sum) :: context
#ifdef BASELINE
  integer, pointer :: fibre(:)
#else
  type(integer_fibre) :: fibre
#endif

  num_vertices = 20000
  num_repetitions = 100
  if (command_argument_count() >= 1) then
     call get_command_argument(1, argument)
     read(argument, *) num_vertices
  end if
  if (command_argument_count() >= 2) then
     call get_command_argument(2, argument)
     read(argument, *) num_repetitions
  end if
  if (num_vertices < degree .or. num_repetitions < 1) error stop 'positive traversal extent required'
  allocate(table(2, degree*num_vertices), tails(degree*num_vertices), heads(degree*num_vertices))
  e = 0
  do i = 1, num_vertices
     do j = 1, degree
        e = e + 1
        tails(e) = i
        heads(e) = modulo(i+j-2, num_vertices) + 1
        table(:,e) = [tails(e), heads(e)]
     end do
  end do
  call domain % declare()
  call sets % bind(domain, counted_set_representation(num_vertices))
  relation = csr_relation('periodic relation', domain, domain, table, sets)
  incidence = stored_directed_graph(num_vertices, tails, heads)
  expected_sum = int(num_repetitions, c_int64_t) * degree * num_vertices * (num_vertices+1_c_int64_t) / 2
  graph_sum = int(num_repetitions, c_int64_t) * (int(e, c_int64_t)*(e+1_c_int64_t)/2 + &
       & degree * int(num_vertices, c_int64_t)*(num_vertices+1_c_int64_t))

  ! A known allocation checks that the link-time counter is active.
  call allocation_begin()
  call relation % image(1, owned)
  num_allocations = allocation_end()
  if (num_allocations < 1 .or. size(owned) /= degree) error stop 'allocation counter did not observe an owning copy'
  print '(a,i0)', 'owning_copy_allocations=', num_allocations

  checksum = 0
  call system_clock(initial_count, count_rate)
  call allocation_begin()
  do repetition = 1, num_repetitions
     do i = 1, num_vertices
#ifdef BASELINE
        fibre => relation % image_view(i)
        do j = 1, size(fibre)
           checksum = checksum + fibre(j)
#else
        fibre = relation % image_view(i)
        do j = 1, fibre % num_members()
           checksum = checksum + fibre % member(j)
#endif
        end do
     end do
  end do
  num_allocations = allocation_end()
  call system_clock(final_count)
  seconds = real(final_count-initial_count, real64) / real(count_rate, real64)
  call report('fibre_index', expected_sum)

  checksum = 0
  call system_clock(initial_count)
  call allocation_begin()
  do repetition = 1, num_repetitions
     do i = 1, num_vertices
#ifdef BASELINE
        fibre => relation % image_view(i)
        call sum_members(fibre)
#else
        fibre = relation % image_view(i)
        call fibre % read(sum_members)
#endif
     end do
  end do
  num_allocations = allocation_end()
  call system_clock(final_count)
  seconds = real(final_count-initial_count, real64) / real(count_rate, real64)
  call report('fibre_read', expected_sum)

  context % value = 0
  call system_clock(initial_count)
  call allocation_begin()
  do repetition = 1, num_repetitions
     do i = 1, num_vertices
#ifdef BASELINE
        fibre => relation % image_view(i)
        call sum_context_members(fibre, context)
#else
        fibre = relation % image_view(i)
        call fibre % read(sum_context_members, context)
#endif
     end do
  end do
  num_allocations = allocation_end()
  call system_clock(final_count)
  seconds = real(final_count-initial_count, real64) / real(count_rate, real64)
  checksum = context % value
  call report('fibre_context', expected_sum)

  checksum = 0
  call system_clock(initial_count)
  call allocation_begin()
  do repetition = 1, num_repetitions
#ifdef BASELINE
     call sum_incoming(incidence % xin, incidence % ein, incidence % tail)
#else
     call incidence % read_incoming(sum_incoming)
#endif
  end do
  num_allocations = allocation_end()
  call system_clock(final_count)
  seconds = real(final_count-initial_count, real64) / real(count_rate, real64)
  call report('incoming_read', graph_sum)

contains

  subroutine report(operation_name, expected)
    character(len=*), intent(in) :: operation_name
    integer(c_int64_t), intent(in) :: expected
    if (checksum /= expected) error stop 'traversal changed incidence values'
    if (num_allocations /= 0) error stop 'traversal allocated storage'
    print '(a,a,f12.6,a,i0,a,i0)', trim(operation_name), ' seconds=', seconds, &
         & ' allocations=', num_allocations, ' checksum=', checksum
  end subroutine report

end program topology_traversal
