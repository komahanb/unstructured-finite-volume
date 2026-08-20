! NAVIGATION PROBE. A function reference as the leftmost part-ref.
! EXPECTED: rejected - 'The leftmost part-ref in a data-ref cannot be
! a function reference'. A function result is not a part-ref, so no
! data-ref may be built on one.
program function_result_base
  use graph_fractal
  implicit none
  integer :: s
  s = null_branch() % status()
  print *, s
end program function_result_base
