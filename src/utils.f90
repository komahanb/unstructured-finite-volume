!=====================================================================!
! Module that contains important vector and matrix operations
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module utils

  use iso_fortran_env, only: dp => real64

  implicit none

  integer, parameter  :: NUM_SPAT_DIM = 3
  real(dp), parameter :: PI = 22.0d0/7.0d0
  real(dp), parameter :: DEG_PER_RAD = 180.0d0/PI ! 180/PI
  real(dp), parameter :: RAD_PER_DEG = PI/180.0d0 ! PI/180

  !-------------------------------------------------------------------!
  ! VECTOR datatype can be used for arrays that are represented in 
  ! three spatial dimensions, for example gravity, acceleration, 
  ! velocity, position, orientation etc
  !-------------------------------------------------------------------!
  
  type vector
     real(dp)  :: x(NUM_SPAT_DIM) = 0.0_dp
   contains
     !procedure :: get
     !procedure :: set
  end type vector

  !-------------------------------------------------------------------!
  ! MATRIX datatype can be used for matrices involved within spatial 
  ! dimensions, for example moment of inertia (J), rotation matrix(C, 
  ! Cdot), angular rates (S, Sdot)
  !-------------------------------------------------------------------!

  type matrix
     real(dp)    :: PSI(NUM_SPAT_DIM, NUM_SPAT_DIM) = 0.0_dp
   contains
     ! procedure :: get
     ! procedure :: set
  end type matrix

  !-------------------------------------------------------------------!
  ! Overload intrinsic norm function to operate on real and complex
  ! numbers
  !-------------------------------------------------------------------!
  interface norm
     module procedure znorm2
     module procedure dnorm2
     module procedure norm_vec
  end interface norm

  !-------------------------------------------------------------------!
  ! Returns the realpart of an array. This overloads the intrinsic
  ! function 'realpart'
  !-------------------------------------------------------------------!
  interface real_part
     module procedure zrealpart
     module procedure drealpart       
  end interface real_part

  !-------------------------------------------------------------------!
  ! Overload * for dot product and other matrix-vector operations
  !-------------------------------------------------------------------!
  ! square(xx), dot(xx, xx), xx*xx       --> all does DOT product
  ! cross product of two vectors a and b --> skew(a)*b
  !
  ! Other scenarios where (*) is used are: 2.0*a, 2.0*A, a*A, A*a, A*A
  !-------------------------------------------------------------------!

  interface operator (*)
     module procedure dot, scal_vec, scal_matrix, &
          & matrix_vector, matrix_matrix, vector_matrix
  end interface operator (*)

  !-------------------------------------------------------------------!
  ! Overload (+) for addition of VECTOR and MATRIX  data types
  !-------------------------------------------------------------------!

  interface operator (+)
     module procedure add_matrices, add_vectors
  end interface operator (+)

  !-------------------------------------------------------------------!
  ! Overload (-) for negation of VECTOR and MATRIX  data types
  !-------------------------------------------------------------------!

  interface operator (-)
     module procedure sub_matrices, sub_vectors, negate_vector, &
          & negate_matrix
  end interface operator (-)

  !-------------------------------------------------------------------!
  ! Gets plain array from VECTOR data type
  !-------------------------------------------------------------------!

  interface array
     module procedure get_array, get_array_1d
  end interface array

  !-------------------------------------------------------------------!
  ! Gets a multidimensional array from MATRIX data type
  !-------------------------------------------------------------------!

  interface matrix
     module procedure new_matrix_from_array, get_matrix, get_matrix_2d
  end interface matrix

contains
  
  !-------------------------------------------------------------------!
  ! Product of a scalar and vector
  !-------------------------------------------------------------------!

  pure elemental function scal_vec(a, b) 

    real(dp), intent (in)      :: a
    type (vector), intent (in) :: b
    type (vector)              :: scal_vec

    scal_vec%x = a*b%x

  end function scal_vec

  !-------------------------------------------------------------------!
  ! Product of a scalar and matrix{aB} = a[B]
  !-------------------------------------------------------------------!

  pure elemental function scal_matrix(a, B) 

    real(dp), intent (in)     :: a
    type(matrix), intent (in) :: B 
    type(matrix)              :: scal_matrix

    scal_matrix%PSI =  a*B%PSI

  end function scal_matrix
  
  !-------------------------------------------------------------------!
  ! Product of a vector and matrix {c} = {a}^T[B]
  !-------------------------------------------------------------------!
  
  pure elemental function vector_matrix(a, B) 

    type(vector), intent (in) :: a
    type(matrix), intent (in) :: B 
    type(vector) ::  vector_matrix

    vector_matrix = vector(matmul(a%x, B%PSI))

  end function vector_matrix

  !-------------------------------------------------------------------!
  ! Product of a matrix with vector {c} = [B]{a}
  !-------------------------------------------------------------------!

  pure elemental function matrix_vector(B, a) 

    type(vector), intent (in) :: a
    type(matrix), intent (in) :: B 
    type(vector) ::  matrix_vector

    matrix_vector = vector(matmul(B % PSI, a % x))

  end function matrix_vector

  !-------------------------------------------------------------------!
  ! Product of two matrices [C] = [A]{B}
  !-------------------------------------------------------------------!

  pure elemental function matrix_matrix(A, B) 

    type(matrix), intent (in) :: A, B 
    type(matrix)              :: matrix_matrix

    matrix_matrix = matrix( matmul(A % PSI, B % PSI))

  end function matrix_matrix

  !-------------------------------------------------------------------!
  ! Returns a skew-symmetric matrix for doing cross product
  !-------------------------------------------------------------------!

  pure elemental function skew(a)

    type(vector), intent(in)           :: a
    type(matrix)                       :: skew
    
    skew = trans( matrix((/ 0.0_dp, a%x(3), -a%x(2), -a%x(3), 0.0_dp, a%x(1), &
         &  a%x(2), -a%x(1), 0.0_dp /)) )

  end function skew

  !-------------------------------------------------------------------!
  ! Returns the vector associated with skew-symmetric matrix 
  !-------------------------------------------------------------------!

  pure elemental function unskew(a)

    type(matrix), intent(in) :: a
    type(vector)             :: unskew

    unskew = vector((/ a%PSI(3,2), a%PSI(1,3), a%PSI(2,1) /))

  end function unskew

  !-------------------------------------------------------------------!
  ! Returns cross product of vectors a and b (can use skew too)
  !-------------------------------------------------------------------!

  pure elemental function cross(a ,b) 

    type(vector), intent(in) :: a, b
    type(vector)             :: cross

    cross = vector(matmul(get_matrix(skew(a)), get_array(b)))

  end function cross

  !-------------------------------------------------------------------!
  ! Returns the dot product of two vectors
  !-------------------------------------------------------------------!

  pure elemental function dot(a,b)

    type (vector), intent (in) :: a, b
    real(dp)                   :: dot

    dot = sum(a%x*b%x) ! a*a

  end function dot

  !-------------------------------------------------------------------!
  ! Returns the magnitude of a vector
  !-------------------------------------------------------------------!

  pure elemental function norm_vec(a)

    type (vector), intent (in) :: a
    real(dp)                   :: norm_vec  

    norm_vec = norm2(a%x) ! dsqrt(dot(a,a)), dsqrt(a*a)

  end function norm_vec

  !-------------------------------------------------------------------!
  ! Get the vector entries as array (can simply use a%x)
  !-------------------------------------------------------------------!

  pure function get_array(a)

    type(vector), intent(in) :: a
    real(dp)                 :: get_array(NUM_SPAT_DIM)

    get_array = a % x

  end function get_array

  !-------------------------------------------------------------------!
  ! Get the entries of a VECTOR of VECTOR as an array
  !-------------------------------------------------------------------!

  pure function get_array_1d(a)

    type(vector), intent(in), dimension(:) :: a

    real(dp)                 :: get_array_1d(NUM_SPAT_DIM*size(a))
    integer                  :: j, is_j, ie_j, n

    n = size(a)
    do j = 1, n
       call split(j,is_j,ie_j) ! split j index storage
       get_array_1d(is_j:ie_j) = get_array(a(j))
    end do

  end function get_array_1d

  !-------------------------------------------------------------------!
  ! Constructor for a new matrix with entries supplied as an array
  ! (only NUM_SPAT_DIM compatible matrices)
  ! -------------------------------------------------------------------!

  pure function new_matrix_from_array(a)

    real(dp), intent(in) :: a(NUM_SPAT_DIM**2)
    type(matrix)         :: new_matrix_from_array

    new_matrix_from_array % PSI = &
         & transpose( reshape(a, (/ NUM_SPAT_DIM, NUM_SPAT_DIM /)) )

  end function new_matrix_from_array

  !-------------------------------------------------------------------!
  ! Get the matrix entries as array
  !-------------------------------------------------------------------!

  pure function get_matrix(A)

    type(matrix), intent(in) :: A
    real(dp)                 :: get_matrix(NUM_SPAT_DIM, NUM_SPAT_DIM)

    get_matrix = A % PSI

  end function get_matrix

  !-------------------------------------------------------------------!
  ! Unwraps a 2d matrix and stores as real numbers
  !-------------------------------------------------------------------!

  pure function get_matrix_2d(A) result(mat)

    type(matrix), dimension(:,:), intent(in) :: A
    real(dp)    :: mat(size(A,1)*NUM_SPAT_DIM, size(A,2)*NUM_SPAT_DIM) 
    integer :: i,j, m, n
    integer :: is_i, ie_i, is_j, ie_j

    n =size(A,1)
    m =size(A,2)

    do i = 1, n
       do j = 1, m
          call split(j,is_j,ie_j) ! split j index storage
          call split(i,is_i,ie_i) ! split i index storage
          mat(is_j:ie_j,is_i:ie_i) = get_matrix(A(j,i))
       end do
    end do

  end function get_matrix_2d


  !-------------------------------------------------------------------!
  ! Unwraps a vector of vector and stores as an array
  !-------------------------------------------------------------------!

  pure function get_vector_elements(a,n)

    integer, intent(in) :: n
    type(vector), intent(in) :: a(n)
    real(dp)    :: get_vector_elements(n*NUM_SPAT_DIM) 
    integer :: i
    integer :: is_i, ie_i

    do i = 1, n
       call split(i,is_i,ie_i) ! split i index storage
       get_vector_elements(is_i:ie_i) = a(i)%x ! get_matrix(A(j,i))
    end do

  end function get_vector_elements

  !-------------------------------------------------------------------!
  ! 
  !-------------------------------------------------------------------!

  pure subroutine split(i,is,ie)

    integer, intent(in) :: i
    integer, intent(out) :: is, ie

    is = (NUM_SPAT_DIM*(i-1) ) + 1
    ie = NUM_SPAT_DIM*( (i-1) + 1)

  end subroutine split

  !-------------------------------------------------------------------!
  ! Transpose of a matrix
  !-------------------------------------------------------------------!

  pure elemental function trans(A)
    
    type(matrix), intent(in) :: A
    type(matrix)             :: trans

    trans = matrix( transpose(A % PSI) )

  end function trans

  !-------------------------------------------------------------------!
  ! Returns the matrix addition of matrices of TYPE matrix
  !-------------------------------------------------------------------!

  pure elemental function add_matrices(A, B)

    type(matrix), intent(in) :: A, B
    type(matrix)  :: add_matrices

    add_matrices%PSI = A%PSI + B%PSI

  end function add_matrices

  !-------------------------------------------------------------------!
  ! Returns the matrix addition of vectors of TYPE matrix
  !-------------------------------------------------------------------!

  pure elemental function add_vectors(a, b)

    type(vector), intent(in) :: a, b
    type(vector)  :: add_vectors

    add_vectors%x = a%x + b%x

  end function add_vectors

  !-------------------------------------------------------------------!
  ! Returns the matrix subtraction of matrices of TYPE matrix
  !-------------------------------------------------------------------!

  pure elemental function sub_matrices(A, B)

    type(matrix), intent(in)  :: A, B
    type(matrix)  :: sub_matrices

    sub_matrices%PSI = A%PSI - B%PSI

  end function sub_matrices

  !-------------------------------------------------------------------!
  ! Returns the matrix addition of vectors of TYPE matrix
  !-------------------------------------------------------------------!

  pure elemental function sub_vectors(a, b)

    type(vector), intent(in) :: a, b
    type(vector)  :: sub_vectors

    sub_vectors%x = a%x - b%x

  end function sub_vectors

  !-------------------------------------------------------------------!
  ! Returns the negative of type matrix
  !-------------------------------------------------------------------!

  pure elemental function negate_matrix(A)

    type(matrix), intent(in)  :: A
    type(matrix)  :: negate_matrix

    negate_matrix%PSI = -A%PSI

  end function negate_matrix

  !-------------------------------------------------------------------!
  ! Returns the negative of type vector
  !-------------------------------------------------------------------!

  pure elemental function negate_vector(a)

    type(vector), intent(in)  :: a
    type(vector)  :: negate_vector

    negate_vector%x = -a%x

  end function negate_vector

  !-------------------------------------------------------------------!
  ! Returns an diagonal matrix of size NUM_SPAT_DIM with the scalar
  ! input as diagonal elements
  ! -------------------------------------------------------------------!

  pure elemental function diag(val)
    
    real(dp), intent(in) :: val
    type(matrix) :: diag
    integer :: k
    
    forall(k=1:NUM_SPAT_DIM)
       diag % PSI (k,k) = val
    end forall

  end function diag

  !-------------------------------------------------------------------!
  ! Generates a nxn identity matrix
  !-------------------------------------------------------------------!

  pure function eye(n)

    integer, intent(in) :: n
    real(dp)    :: eye(n,n)
    integer     :: i, j

    ! generate a zero matrix
    eye = 0.0_dp
    
    ! replace the diagonals with one
    forall(i=1:n)
       eye(i,i) = 1.0_dp
    end forall

  end function eye

  !-------------------------------------------------------------------!
  ! Convert from degree to radian
  !-------------------------------------------------------------------!

  pure elemental function deg2rad(deg)

    real(dp), intent(in) :: deg
    real(dp) :: deg2rad

    deg2rad  =  deg*RAD_PER_DEG

  end function deg2rad

  !-------------------------------------------------------------------!
  ! Convert from radian to degree
  !-------------------------------------------------------------------!

  pure elemental function rad2deg(rad)

    real(dp), intent(in)  :: rad
    real(dp) :: rad2deg

    rad2deg  =  rad*DEG_PER_RAD

  end function rad2deg

  !===================================================================!
  ! Norm of a complex number array
  !===================================================================!

  real(dp) pure function znorm2(z)

    complex(dp), dimension(:), intent(in) :: z
    integer :: j, n

    znorm2 = 0.0d0
    n = size(z)
    do j = 1, n
       znorm2 = znorm2 + dsqrt(dble(z(j))**2 + aimag(z(j))**2)
    end do

  end function znorm2

  !===================================================================!
  ! Norm of a complex real number array
  !===================================================================!
  
  real(dp) pure function dnorm2(z) result(val)
    
    real(dp), dimension(:), intent(in) :: z

    val = norm2(z)

  end function dnorm2


  !===================================================================!
  ! Norm of a complex number array
  !===================================================================!

  real(dp) pure elemental function zrealpart(z) result(val) 
    
    complex(dp), intent(in) :: z
    
    ! Use the intrinsic function to get the realpart
    val = real(z)
    
  end function zrealpart

  !===================================================================!
  ! Norm of a complex real number array
  !===================================================================!
  
  real(dp) pure elemental function drealpart(z) result(val)
    
    real(dp), intent(in) :: z

    ! Return the real number as it is
    val = z
    
  end function drealpart

end module utils
