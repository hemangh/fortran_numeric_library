module dvr
implicit NONE

!integer, parameter :: dp = selected_real_kind(16, 99)


type dvr_type
   !number of polynomials
   integer :: num
   !points of grid
   real(8),allocatable :: pt(:)
   !weights of polynomials
   real(8),allocatable :: w(:)
   !Legendre polynomials on the grid points
   real(8),allocatable :: lp(:,:)
   !first derivative of ploynomials 
   real(8), allocatable :: dlp(:,:)
   !second derivative of polynomials
   real(8), allocatable :: ddlp(:,:)
   !if either end is fixed but then dropped (zeroed) in calculations
   logical::   dropped(2)
   contains
     procedure :: initilize_dvr_type
     generic   :: initialize => initilize_dvr_type 
end type

type cont_dvr_type
   !number of (radial) grid points
   integer                    :: num
   !normalized polynomial at grid points
   real(8),allocatable       :: dvr(:), ddvr(:)
   ! constituent dvr_types
   type(dvr_type),allocatable :: dvr_elem(:)
   integer                    :: num_regions
   contains
    procedure  :: initialize_cont_dvr_type
    generic    :: initialize => initialize_cont_dvr_type 
end type  

contains

  subroutine initilize_dvr_type(self,grid_pts,grid_weights,dropped) 
    class(dvr_type) :: self
    real(8)        :: grid_pts(:), grid_weights(:)
    integer         :: gsize
    logical, optional         :: dropped(2)
    !allocate
    gsize    = size(grid_pts)
    self%num = gsize
    if(self%num < 3) stop 'cannot make dvrs with less than 3 grid points'
    
    allocate(self%pt(gsize),self%w(gsize))
    allocate(self%lp(gsize,gsize),self%dlp(gsize,gsize),self%ddlp(gsize,gsize))
    
    self%pt = grid_pts
    self%w  = grid_weights

    call cpoly(self%lp, self%dlp, self%ddlp, grid_pts, gsize)

    if (present(dropped)) then
      self%dropped = dropped
    end if  
         
  end subroutine
  
  subroutine initialize_cont_dvr_type(self,dvr_arr) 
    class(cont_dvr_type) :: self
    type(dvr_type)       :: dvr_arr(:)
    end subroutine
!

! SUBROUTINE lgngr(p,dp,ddp,x,y,nx,ny,type,drctv)
! !***begin prologue     lgngr
! !***date written       940504   (yymmdd)
! !***revision date               (yymmdd)
! !***keywords
! !***author             schneider, barry (nsf)
! !***source             %W% %G%
! !***purpose            lagrange polynomials at arbitrary points.
! !***description
! !***
! !
! !
! !***references
! !
! !***routines called
! !
! !***end prologue       lgngr
! !
! IMPLICIT NONE
! REAL*8, DIMENSION (ny,nx)           :: p
! REAL*8, DIMENSION (ny,nx)           :: dp
! REAL*8, DIMENSION (ny,nx)           :: ddp
! REAL*8, DIMENSION(nx)               :: x
! REAL*8, DIMENSION(ny)               :: y
! REAL*8, DIMENSION(:), ALLOCATABLE   :: xt
! REAL*8, DIMENSION(:), ALLOCATABLE   :: yt
! REAL*8                              :: sn
! REAL*8                              :: ssn
! REAL*8                              :: fac
! CHARACTER (LEN = *), OPTIONAL       :: drctv
! CHARACTER (LEN = *), OPTIONAL       :: type
! INTEGER                             :: nx
! INTEGER                             :: ny
! INTEGER                             :: i
! INTEGER                             :: j
! INTEGER                             :: k
! INTEGER                             :: first
! INTEGER                             :: second
! INTEGER                             :: zerfac

! !
! !     generate polynomials and derivatives with respect to x
! !
! p(:,:) = 1.d0
! IF (present(type) ) THEN
!    ALLOCATE(xt(nx),yt(ny))
!    xt(:) = x(:)
!    yt(:) = y(:)
!    x(:) = x(:) * x(:)
!    y(:) = y(:) * y(:)
! END IF
! DO i=1,ny
!   zerfac = 0
!   DO j=1,nx
!      fac =  y(i) - x(j)
!      IF(abs(fac) <= 1.d-10) THEN
!         zerfac = j
!      ENDIF
!   END DO
!   DO j=1,nx
!      DO k = 1, j-1
!         p(i,j) = p(i,j) * ( y(i) - x(k) )   &
!                         / ( x(j) - x(k) )
!      END DO
!      DO k=j+1,nx
!         p(i,j) = p(i,j) * ( y(i) - x(k) )   &
!                         / ( x(j) - x(k) )
!      END DO
!      IF(present(drctv) ) THEN
!          IF ( abs(p(i,j)) > 1.d-10) THEN
!               sn = 0.d0
!               ssn = 0.d0
!               DO k=1,j-1
!                 !if(i /= k) then
!                  fac = 1.d0/( y(i) - x(k) )
!                  sn = sn + fac
!                  ssn = ssn + fac*fac
!                 !endif
!               END DO
!               DO k=j+1,nx
!                 !if(i /= k) then
!                  fac = 1.d0/( y(i) - x(k) )
!                  sn = sn + fac
!                  ssn = ssn + fac*fac
!                 !end if
!               END DO
!               dp(i,j) = sn*p(i,j)
!               ddp(i,j) = sn*dp(i,j) - ssn*p(i,j)
!          ELSE
!               first=j
!               second=zerfac
!               IF(first > second) THEN
!                  first=zerfac
!                  second=j
!               END IF
!               sn = 1.d0
!               ssn = 0.d0
!               DO k=1,first-1
!                  fac = 1.d0/( x(j) - x(k) )
!                  sn = sn*fac*( y(i) - x(k) )
!                  ssn = ssn + 1.d0/(y(i) - x(k))
!               END DO
!               DO k=first+1,second-1
!                  fac = 1.d0/( x(j) - x(k) )
!                  sn = sn*fac*( y(i) - x(k) )
!                  ssn = ssn + 1.d0/( y(i) - x(k) )
!               END DO
!               DO k=second+1,nx
!                  fac = 1.d0/( x(j) - x(k) )
!                  sn = sn*fac*( y(i) - x(k) )
!                  ssn = ssn + 1.d0/( y(i) - x(k) )
!               END DO
!               dp(i,j) = sn/( x(j) - x(zerfac) )
!               ddp(i,j) = 2.d0*ssn*dp(i,j)
!          END IF
!      END IF
! !
!   END DO
! END DO
! !
! IF (present(type)) THEN
!    DO i=1,ny
!       ddp(i,:) = 2.d0*dp(i,:) + 4.d0 * yt(i) * yt(i) * ddp(i,:)
!       dp(i,:) = 2.d0 * yt(i) * dp(i,:)
!    END DO
!    x(:) = xt(:)
!    y(:) = yt(:)
!    DEALLOCATE(xt,yt)
! !
! END IF

! END SUBROUTINE Lgngr

! SUBROUTINE cpoly(cp,dcp,ddcp,pt,n)
! IMPLICIT NONE
! INTEGER                             :: n
! REAL(8), DIMENSION(:,:)           :: cp
! REAL(8), DIMENSION(:,:)           :: dcp
! REAL(8), DIMENSION(:,:)           :: ddcp
! REAL(8), DIMENSION(:)             :: pt

! call lgngr(cp,dcp,ddcp,pt,pt,n,n,drctv='on')

! end subroutine cpoly


! !*******************************************************************************
! !> map [-1,1] quadrature points on arbitrary [a,b]=[endpts(1), endpts(2)]
! !>  y = (b-a)/2 x + (a+b)/2, with x the [-1,1] region quadrature points and
! !>  y the [a,b] region quadrature points and
! !>  w_[a,b] = (b-a)* w_[-1,1]/2 are the mapped weights
! !*******************************************************************************
!   subroutine map_to_endpoints( t, w, endpts, n )
!   !>all rules are from -1 to 1 range and have to be mapped to
!   !>the desired range:

!     implicit none
!     integer                                       :: i
!     integer                                       :: n
!     real(8), dimension(n)                         :: t
!     real(8), dimension(n)                         :: w
!     real(8), dimension(2)                         :: endpts

!   do i= 1 , n
!      t(i) = t(i) * (endpts(2)-endpts(1)) * 0.5 + (endpts(2)+endpts(1)) * 0.5
!      w(i) = w(i) * 0.5 * ( endpts(2) - endpts(1) )
!   end do
!   end subroutine map_to_endpoints

! !*******************************************************************************
! subroutine function_class(kind, n, b, a, muzero, alpha, beta)
! !*******************************************************************************
! !>           this procedure supplies the coefficients a(j), b(j) of the
! !>        recurrence relation
! !>
! !>             b p (x) = (x - a ) p   (x) - b   p   (x)
! !>              j j            j   j-1       j-1 j-2
! !>
! !>        for the various classical (normalized) orthogonal polynomials,
! !>        and the zero-th moment
! !>
! !>             muzero = integral w(x) dx
! !>
! !>        of the given polynomial   weight function w(x).  since the
! !>        polynomials are orthonormalized, the tridiagonal matrix is
! !>        guaranteed to be symmetric.
! !>
! !>           the input parameter alpha is used only for laguerre and
! !>        jacobi polynomials, and the parameter beta is used only for
! !>        jacobi polynomials.  the laguerre and jacobi polynomials
! !>        require the gamma function.
! !*******************************************************************************
!     character(len=*)           :: kind
!     integer,intent(in)         :: n
!     integer                    :: nm1, i
!     real(8)                    :: muzero
!     real(8),optional           :: alpha, beta
!     real(8),dimension(n)       :: a,b
!     real(8),parameter          :: pi = atan(1.d0)*4.d0
!     real(8)                    :: a2b2,ab, abi
! !    real(8)                    :: gamfun use instead the standard gamma(x) of compiler
!     nm1 = n - 1

!     select case (kind)

!         case('legendre')
!         !>      legendre polynomials p(x)
!         !>         on (-1, +1), w(x) = 1.
!            muzero = 2.0d0
!            do i = 1, nm1
!               a(i) = 0.0d0
!               abi = real(i)
!               b(i) = abi/ sqrt(4.d0*abi*abi - 1.0d0  )
!            end do
!            a(n) = 0.0d0

!         case('chebyshev1')
!              !>       chebyshev polynomials of the first kind t(x)
!              !>        on (-1, +1), w(x) = 1 / sqrt(1 - x*x)
!              !>
!                muzero = pi
!                do i = 1, nm1
!                     a(i) = 0.0d0
!                     b(i) = 0.5d0
!               end do
!               b(1) =  sqrt(0.5d0  )
!               a(n) = 0.0d0

!         case('chebyshev2')
!                 !>    chebyshev polynomials of the second kind u(x)
!                 !>     on (-1, +1), w(x) = sqrt(1 - x*x)
!                 !>
!                 muzero = pi/2.0d0
!                 do i = 1, nm1
!                     a(i) = 0.0d0
!                     b(i) = 0.5d0
!                 end do
!                 a(n) = 0.0d0

!          case('hermite')
!                 !>       hermite polynomials h(x) on (-infinity,
!                 !>       +infinity), w(x) = exp(-x**2)
!                 muzero =  sqrt(pi)
!                 do i = 1, nm1
!                     a(i) = 0.0d0
!                     b(i) =  sqrt(i/2.0d0  )
!                 end do
!                 a(n) = 0.0d0

! !next two cases require a gamma function definition as well as alpha beta
!                   case('jacobi')
!                     !>     jacobi polynomials p(alpha, beta)(x) on
!                     !>   (-1, +1), w(x) = (1-x)**alpha + (1+x)**beta, alpha and
!                     !>        beta greater than -1

!                     if (present(alpha).and.present(beta))then

!                       ab = alpha + beta
!                       abi = 2.0d0   + ab
!                       ! muzero = 2.0d0   ** (ab + 1.0d0  ) * gamfun(alpha + 1.0d0  ) * &
!                       ! gamfun(beta + 1.0d0  ) / gamfun(abi)
!                       muzero = 2.0d0   ** (ab + 1.0d0  ) * gamma(alpha + 1.0d0  ) * &
!                       gamma(beta + 1.0d0  ) / gamma(abi)                      
!                       a(1) = (beta - alpha)/abi
!                       b(1) =  sqrt(4.0d0  *(1.0d0  + alpha)*(1.0d0   + beta)/ &
!                              ((abi + 1.0d0  )*  abi*abi))
!                       a2b2 = beta*beta - alpha*alpha
!                       do i = 2, nm1
!                         abi = 2.0d0  *i + ab
!                         a(i) = a2b2/((abi - 2.0d0  )*abi)
!                         b(i) =  sqrt (4.0d0  *i*(i + alpha)*(i + beta)*(i + ab)/ &
!                                  ((abi*abi - 1)*abi*abi))
!                       end do
!                       abi = 2.0d0  *n + ab
!                       a(n) = a2b2/((abi - 2.0d0  )*abi)

!                     else

!                       stop "jacobi functions need alpha and beta optionals"

!                     end if




!                   case('lagurre')
!                     !>     laguerre polynomials l(alpha)(x) on
!                     !>  (0, +infinity), w(x) = exp(-x) * x**alpha, alpha greater
!                     !>   than -1.
!                     !>

!                     if (present(alpha))then

!                       ! muzero = gamfun(alpha + 1.0d0  )
!                       muzero = gamma(alpha + 1.0d0  )
!                       do  i = 1, nm1
!                              a(i) = 2.0d0  *i - 1.0d0   + alpha
!                              b(i) =  sqrt(i*(i + alpha))
!                       end do
!                       a(n) = 2.0d0  * n - 1 + alpha


!                   ! case('simpson')

!                   ! case('trapezoidal')

!                     else

!                       stop "lagurre functions need alpha optional"

!                     endif

!          case default
!                  stop 'wrong name for function_class in module dvr'

!    end select

! end subroutine function_class

! !******************************************************************************
! subroutine gauss_quadrature(kind, n, kpts, endpts, t, w, alpha, beta)
!   !>        n        the number of points used for the quadrature rule
!   !>        alpha    real parameter used only for gauss-jacobi and gauss-
!   !>                 laguerre quadrature optional arguments.
!   !>        beta     real parameter used only for gauss-jacobi quadrature--
!   !>                 optional arguments.
!   !>        kpts     (integer) normally 0, unless the left or right end-
!   !>                 point (or both) of the interval is required to be a
!   !>                 node (this is called gauss-radau or gauss-lobatto
!   !>                 quadrature).  then kpts is the number of fixed
!   !>                 endpoints (1 or 2).
!   !>        endpts   real array of length 2.  contains the values of
!   !>                 any fixed endpoints, if kpts = 1 or 2.
!   !>        b        real scratch array of length n
!   !>
!   !>     output parameters (both arrays of length n)
!   !>
!   !>        t        will contain the desired nodes x(1),,,x(n)
!   !>        w        will contain the desired weights c(1),,,c(n)
! !******************************************************************************

!       integer, intent(in)              :: n
!       integer, intent(in)              :: kpts
!       integer                          :: ierr, i, nm1
!       real(8),dimension(2)             :: endpts
!       real(8),dimension(n)             :: b
!       real(8),dimension(n),intent(out) :: t,w
!       real(8)                          :: muzero,gam,t1
!       real(8)                          :: h
!       real(8)                          :: gbslve
!       real(8), optional                :: alpha, beta
!       character(len=*)                 :: kind
!       !print*, trim(kind)

!       if (trim(kind) == 'simpson' .or. trim(kind)=='trapezoidal')then

!        !case of simpson or trapezoidal kpts does not matter here.
!         ! The two ends are always fixed or kpts=2
!          select case(trim(kind))

!          case('simpson')
!            if(mod(n,2)==0) then
!              !c         call lnkerr('n must be odd for simpson rule')
!              stop 'n must be odd for simpson rule'
!              !n = n + 1
!            endif
!            if(n.le.1) then
!              t(1) = 0.d+00
!              w(1) = 2.d+00
!            endif
!            h = 2.d+00/(n-1)
!            t(1) = -1.d+00
!            t(n) = 1.d+00
!            w(1) = h/3.d+00
!            w(n) = h/3.d+00
!            nm1 = n-1
!            do i=2,nm1
!              t(i) = t(i-1) + h
!              w(i) = 4.d+00 - 2.d+00*(i-2*(i/2))
!              w(i) = w(i)*h/3.d+00
!            end do

!          case('trapezoidal')

!            h   = 2.d+00/(n-1)
!            w(1)= h*0.5d+00
!            w(n)= h*0.5d+00
!            t(1)= -1.d+00
!            t(n)=  1.d+00
!            do i=2,n-1
!               w(i)= h
!               t(i)= t(i-1) + h
!            end do

!          case default

!            write(*,*) 'kind:',kind
!            stop 'quadrature type is not defined'

!          end select

!       else

!          if (present(alpha).and.present(beta))then

!             call function_class (trim(kind), n, b, t, muzero,alpha, beta)

!          elseif(present(alpha))then

!             call function_class (trim(kind), n, b, t, muzero,alpha)

!          else

!             call function_class (trim(kind), n, b, t, muzero)

!          endif


! !
! !>           the matrix of coefficients is assumed to be symmetric.
! !>          the array t contains the diagonal elements, the array
! !>           b the off-diagonal elements.
! !>          make appropriate changes in the lower right 2 by 2
! !>           submatrix.
! !

!         if(kpts.eq.1)then

!         !
!         !>           if kpts=1, only t(n) must be changed
!         !
!               t(n) =gbslve(endpts(1), n, t, b)*b(n-1)**2 + endpts(1)


!         elseif(kpts.eq.2)then

!         !
!         !>           if kpts=2, t(n) and b(n-1) must be recomputed
!         !
!           gam =gbslve(endpts(1), n, t, b)
!           t1 = ((endpts(1) - endpts(2))/(gbslve(endpts(2), n, t, b) - gam))
!           b(n-1) =  sqrt(t1)
!           t(n) = endpts(1) + gam*t1


!         elseif (kpts .ne. 0) then

!            stop 'kpts in gauss_quadrature routine in dvr module has invalid value'

!         end if

! !
! !           note that the indices of the elements of b run from 1 to n-1
! !           and thus the value of b(n) is arbitrary.
! !           now compute the eigenvalues of the symmetric tridiagonal
! !           matrix, which has been modified as necessary.
! !           the method used is a ql-type method with origin shifting
! !
!         w(1) = 1.0d0
!         do  i = 2, n
!            w(i) = 0.0d0
!         end do
! !
!         call gbtql2 (n, t, b, w, ierr)

!         do  i = 1, n
!            w(i) = muzero * w(i) * w(i)
!         end do


!     end if

! !

! end subroutine gauss_quadrature


end module dvr
