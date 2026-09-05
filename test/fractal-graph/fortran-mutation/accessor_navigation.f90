! CANDIDATES C and D, navigation. With branch private, g % branch(i)
! is a function reference, so g % branch(i) % status() chains past one.
! EXPECTED: rejected. gfortran words the diagnostic by form - here
! "Invalid character in name", in chained_navigation.f90
! "Unclassifiable statement" - but the rule is one: a function result
! is not a part-ref. The mechanism that makes C immutable is therefore
! the same rule that removes direct navigation; they are not separable.
program accessor_navigation
  use accessor_candidates
  implicit none
  type(graph_t) :: g
  integer       :: s
  s = g % by_value(1) % status()
  print *, s
end program accessor_navigation
