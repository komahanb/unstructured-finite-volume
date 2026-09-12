!=====================================================================!
! Listed subsets: construction, membership and field transport.
!
! A carrier of N members, a subset of M members chosen by a fixed
! pseudo-random permutation, declared ascending, descending or in
! permutation order, and optionally with every member repeated once.
!
!      construct   listed_set_representation(values), membership by
!                  local_index for members and for non-members
!      transport   a field on the subset of a chain graph's vertices,
!                  restricted to P linear parts and assembled back
!
! The oracle is the marker-array construction: first occurrences in
! declaration order, and membership by a logical array over 1..N.
!
! usage: subset mode N M order duplicated parts repetitions
!=====================================================================!

program subset_scaling

  use iso_fortran_env       , only : int64
  use util_precision        , only : dp
  use graph_fractal         , only : graph
  use map_set_representation, only : counted_set_representation, listed_set_representation
  use map_set_store         , only : set_store
  use view_directed         , only : directed_graph
  use view_directed_stored  , only : stored_directed_graph
  use field_calculus        , only : field
  use field_stored          , only : stored_field
  use transform_partitioner , only : partitioner, PARTITION_LINEAR
  use transform_assembler   , only : assembler
  use relation_partition    , only : partition_relation
  use benchmark_measurement , only : phase, begin_phase, end_phase, record, peak_rss_kilobytes, &
       &                             argument_integer, argument_string, pseudo_random_permutation

  implicit none

  character(len=:), allocatable :: mode, order, tokens
  character(len=128) :: written
  integer :: num_carrier, num_subset, parts, repetitions, duplicated, i, k, num_failures, total
  integer(int64) :: position_sum
  integer, allocatable :: permutation(:), members(:), input(:), expected(:), complement(:), positions(:)
  logical, allocatable :: member(:), assembled(:)
  type(listed_set_representation) :: representation
  type(phase) :: measured
  logical :: satisfied

  mode        = argument_string(1, 'construct')
  num_carrier = argument_integer(2, 100000)
  num_subset  = argument_integer(3, 10000)
  order       = argument_string(4, 'random')
  duplicated  = argument_integer(5, 0)
  parts       = argument_integer(6, 4)
  repetitions = argument_integer(7, 1)
  if (num_subset > num_carrier .or. num_subset < 0 .or. num_carrier < 1) then
     error stop 'subset: the subset lies within the carrier'
  end if
  write(written, '(a,a,a,i0,a,i0,a,a,a,i0,a,i0)') 'suite=subset mode=', mode, ' carrier=', num_carrier, &
       & ' subset=', num_subset, ' order=', order, ' duplicated=', duplicated, ' parts=', parts
  tokens = trim(written)
  num_failures = 0

  permutation = pseudo_random_permutation(num_carrier, 7919)
  allocate(member(num_carrier), source=.false.)
  member(permutation(1:num_subset)) = .true.
  select case (order)
  case ('ascending')
     members = pack([(i, i = 1, num_carrier)], member)
  case ('descending')
     members = pack([(i, i = 1, num_carrier)], member)
     members = members(size(members):1:-1)
  case ('random')
     members = permutation(1:num_subset)
  case default
     error stop 'subset: the order is ascending, descending or random'
  end select
  if (duplicated /= 0) then
     input = [members, members]
  else
     input = members
  end if
  ! the oracle: first occurrences in declaration order
  allocate(positions(num_carrier), source=0)
  allocate(expected(size(input)))
  total = 0
  do i = 1, size(input)
     if (positions(input(i)) /= 0) cycle
     total = total + 1
     positions(input(i)) = total
     expected(total) = input(i)
  end do
  expected = expected(1:total)
  complement = pack([(i, i = 1, num_carrier)], .not. member)
  if (size(complement) > num_subset) complement = complement(1:num_subset)
  if (size(complement) == 0) complement = [(num_carrier + i, i = 1, num_subset)]

  select case (mode)
  case ('construct')
     call begin_phase(measured)
     do k = 1, repetitions
        representation = listed_set_representation(input)
     end do
     call end_phase(measured)
     call record(tokens, 'construction', measured)

     satisfied = representation % num_members() == total
     if (satisfied) then
        do i = 1, total
           satisfied = satisfied .and. representation % member(i) == expected(i)
        end do
        do i = 1, num_carrier
           satisfied = satisfied .and. representation % local_index(i) == positions(i)
        end do
     end if
     call verified(satisfied, 'members_and_inverse_equal_oracle')

     position_sum = 0
     call begin_phase(measured)
     do k = 1, repetitions
        do i = 1, size(members)
           position_sum = position_sum + representation % local_index(members(i))
        end do
     end do
     call end_phase(measured)
     call record(tokens, 'membership_members', measured)
     call verified(position_sum == int(repetitions, int64) * (int(total, int64) * int(total + 1, int64)) / 2, &
          & 'member_positions_sum')

     position_sum = 0
     call begin_phase(measured)
     do k = 1, repetitions
        do i = 1, size(complement)
           position_sum = position_sum + representation % local_index(complement(i))
        end do
     end do
     call end_phase(measured)
     call record(tokens, 'membership_non_members', measured)
     if (position_sum /= 0) then
        do i = 1, size(complement)
           if (representation % local_index(complement(i)) /= 0) then
              write(*, '(a,i0,a,i0,a,i0)') 'non-member ', complement(i), ' has position ', &
                   & representation % local_index(complement(i)), ' and oracle position ', positions(complement(i))
              exit
           end if
        end do
     end if
     call verified(position_sum == 0, 'non_members_have_no_position')

  case ('transport')
     call transported()

  case default
     error stop 'subset: the mode is construct or transport'
  end select

  write(*, '(a,1x,a,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0)') 'summary', tokens, 'input=', size(input), &
       & 'retained=', total, 'queries=', size(complement), 'peak_rss_kilobytes=', peak_rss_kilobytes()
  if (num_failures > 0) error stop 'subset: a verification failed'

contains

  subroutine verified(condition, name)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: name
    if (.not. condition) num_failures = num_failures + 1
    write(*, '(a,1x,a,1x,a,a,1x,a,l1)') 'verification', tokens, 'name=', name, 'satisfied=', condition
  end subroutine verified

  ! A chain graph on the carrier; the subset field restricted to every
  ! linear part and assembled back onto a subobject of the whole carrier.
  subroutine transported()
    type(stored_directed_graph) :: g
    type(set_store) :: sets
    type(graph) :: carrier, chosen, transported_domain
    type(stored_field) :: source_field
    type(partitioner) :: p
    type(assembler) :: a
    type(partition_relation) :: rel
    class(directed_graph), allocatable :: part
    class(field), allocatable :: part_data, assembled_data
    real(dp), allocatable :: values(:)
    integer, allocatable :: mapped(:)
    type(phase) :: cut, restriction, assembly
    integer :: v

    g = stored_directed_graph(num_carrier, tails=[(v, v = 1, num_carrier - 1)], heads=[(v + 1, v = 1, num_carrier - 1)])
    carrier = g % vertex_set()
    call sets % bind(carrier, counted_set_representation(num_carrier))
    call begin_phase(measured)
    call sets % declare_subobject(chosen, input, 'chosen', carrier)
    call end_phase(measured)
    call record(tokens, 'declaration', measured)
    source_field = stored_field('q', chosen, total)
    call source_field % set_real_vector([(10.0_dp * expected(i) + 1.0_dp, i = 1, total)])

    allocate(assembled(num_carrier), source=.false.)
    satisfied = .true.
    cut % seconds = 0.0_dp
    restriction % seconds = 0.0_dp
    assembly % seconds = 0.0_dp
    do k = 1, parts
       call begin_phase(measured)
       p = partitioner(PARTITION_LINEAR, num_parts=parts, part=k)
       call p % partition_graph(g, part, rel)
       call end_phase(measured)
       call accumulate(cut, measured)

       call begin_phase(measured)
       call sets % bind(part % vertex_set(), counted_set_representation(part % num_vertices()))
       call sets % bind(part % edge_set(), counted_set_representation(part % num_edges()))
       call p % partition_data(rel, g, source_field, part, sets, part_data)
       call end_phase(measured)
       call accumulate(restriction, measured)

       call begin_phase(measured)
       call a % assemble_data(rel, part, part_data, g, sets, assembled_data)
       call end_phase(measured)
       call accumulate(assembly, measured)

       transported_domain = assembled_data % domain()
       call sets % members_of(transported_domain, mapped)
       call assembled_data % real_vector(values)
       satisfied = satisfied .and. size(values) == size(mapped)
       do i = 1, size(mapped)
          v = mapped(i)
          if (v < 1 .or. v > num_carrier) then
             satisfied = .false.
             cycle
          end if
          satisfied = satisfied .and. member(v) .and. .not. assembled(v)
          assembled(v) = .true.
          if (i <= size(values)) satisfied = satisfied .and. values(i) == 10.0_dp * v + 1.0_dp
       end do
    end do
    satisfied = satisfied .and. all(assembled .eqv. member)
    call record(tokens, 'partition_graph', cut)
    call record(tokens, 'partition_data', restriction)
    call record(tokens, 'assemble_data', assembly)
    call verified(satisfied, 'assembled_members_and_values_equal_source')
  end subroutine transported

  subroutine accumulate(total_phase, one)
    type(phase), intent(inout) :: total_phase
    type(phase), intent(in) :: one
    total_phase % seconds = total_phase % seconds + one % seconds
    total_phase % calls = total_phase % calls + one % calls
    total_phase % bytes = total_phase % bytes + one % bytes
  end subroutine accumulate

end program subset_scaling
