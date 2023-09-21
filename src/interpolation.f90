module interpolation

implicit none
private
public spline, ispline, &
       spline_type, linear_intp!intp_constants, spline_interpolate, line_interpolate


type spline_type
     integer  :: dat_size
     integer  :: order
     real(8), allocatable :: x(:)
     real(8), allocatable :: y(:)
     real(8), allocatable :: const(:,:)
     logical               :: called = .false.
   contains
     procedure :: initialize
     procedure  :: evaluate
end type

contains

subroutine initialize(sp, x, y)
    class(spline_type) :: sp
    real(8)          :: x(:), y(:)
    integer           :: order
    integer           :: nsize

    nsize = size(x)
    sp%dat_size = nsize

    if(nsize <=1 .or. nsize /= size(y))then
      stop 'wrong data size or mismatched arrays'
      print*, 'x size:', nsize
      print*, 'y size:', size(y)
    end if

    order = 3 ! always cubic spline
    sp%order = order

    if (sp%called .eqv. .true.)then
      deallocate(sp%x, sp%y, sp%const)
    end if

    allocate(sp%x(nsize), sp%y(nsize), sp%const(nsize,order))

    sp%x = x
    sp%y = y

    call spline (x, y, sp%const(:,1), sp%const(:,2), sp%const(:,3), nsize)

    sp%called = .true.

end subroutine


function evaluate(sp,x) result(y)
     class(spline_type) :: sp
     real(8) :: x
     real(8) :: y

     if (sp%called .eqv. .false.) then
       stop 'un-initialized spline interpolation is asked for'
     end if

      y = ispline(x, sp%x, sp%y, sp%const(:,1), sp%const(:,2), &
                                   sp%const(:,3), sp%dat_size)

end function
!**************************************************************************************
   subroutine spline (x, y, b, c, d, n)
!======================================================================
!  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
!  for cubic spline interpolation
!  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!  for  x(i) <= x <= x(i+1)
!  Alex G: January 2010
!----------------------------------------------------------------------
!  input..
!  x = the arrays of data abscissas (in strictly increasing order)
!  y = the arrays of data ordinates
!  n = size of the arrays xi() and yi() (n>=2)
!  output..
!  b, c, d  = arrays of spline coefficients
!  comments ...
!  spline.f90 program is based on fortran version of program spline.f
!  the accompanying function fspline can be used for interpolation
!======================================================================

integer n
double precision x(n), y(n), b(n), c(n), d(n)
integer i, j, gap
double precision h

gap = n-1
! check input
if ( n < 2 ) return
if ( n < 3 ) then
  b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
  c(1) = 0.
  d(1) = 0.
  b(2) = b(1)
  c(2) = 0.
  d(2) = 0.
  return
end if
!
! step 1: preparation
!
d(1) = x(2) - x(1)
c(2) = (y(2) - y(1))/d(1)
do i = 2, gap
  d(i) = x(i+1) - x(i)
  b(i) = 2.0*(d(i-1) + d(i))
  c(i+1) = (y(i+1) - y(i))/d(i)
  c(i) = c(i+1) - c(i)
end do
!
! step 2: end conditions
!
b(1) = -d(1)
b(n) = -d(n-1)
c(1) = 0.0
c(n) = 0.0
if(n /= 3) then
  c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
  c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
  c(1) = c(1)*d(1)**2/(x(4)-x(1))
  c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
end if
!
! step 3: forward elimination
!
do i = 2, n
  h = d(i-1)/b(i-1)
  b(i) = b(i) - h*d(i-1)
  c(i) = c(i) - h*c(i-1)
end do
!
! step 4: back substitution
!
c(n) = c(n)/b(n)
do j = 1, gap
  i = n-j
  c(i) = (c(i) - d(i)*c(i+1))/b(i)
end do
!
! step 5: compute spline coefficients
!
b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
do i = 1, gap
  b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
  d(i) = (c(i+1) - c(i))/d(i)
  c(i) = 3.*c(i)
end do
c(n) = 3.0*c(n)
d(n) = d(n-1)
end subroutine spline

  function ispline(u, x, y, b, c, d, n)
!======================================================================
! function ispline evaluates the cubic spline interpolation at point z
! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
! where  x(i) <= u <= x(i+1)
!----------------------------------------------------------------------
! input..
! u       = the abscissa at which the spline is to be evaluated
! x, y    = the arrays of given data points
! b, c, d = arrays of spline coefficients computed by spline
! n       = the number of data points
! output:
! ispline = interpolated value at point u
!=======================================================================

double precision ispline
integer n
double precision  u, x(n), y(n), b(n), c(n), d(n)
integer i, j, k
double precision dx

! if u is ouside the x() interval take a boundary value (left or right)
if(u <= x(1)) then
  ispline = y(1)
  return
end if
if(u >= x(n)) then
  ispline = y(n)
  return
end if

!*
!  binary search for for i, such that x(i) <= u <= x(i+1)
!*
i = 1
j = n+1
do while (j > i+1)
  k = (i+j)/2
  if(u < x(k)) then
    j=k
    else
    i=k
   end if
end do
!*
!  evaluate spline interpolation
!*
dx = u - x(i)
ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
end function ispline
!!******************************************************
function linear_intp(u,x,y,n)result(uy)
implicit none
integer n
double precision x(n), y(n), u, uy, slope
integer i, j, flag

!determine the interval
!!print*,u

flag= 0
do i=1,n-1
    if(u == x(i))then
       uy=y(i)
       flag= 1
       exit
    end if
if(x(i) .lt. u .AND. u .lt. x(i+1))exit
!!	print*,i, x(i), x(i+1)
end do
!check last point
if(u==x(n))then
  uy=y(n)
  flag=1
end if
!!print*,x(i),x(i+1), y(i), y(i+1)
if (flag .ne. 1) then
!!    print*,'i=',i
if(y(i+1)-y(i)/=0.d0)then
    slope=(y(i+1)-y(i))/(x(i+1)-x(i))
    uy=y(i)+slope*(u-x(i))
  else
    uy=y(i)
  end if
end if
!!print*,uy

end function
!*********************************************
end module interpolation
