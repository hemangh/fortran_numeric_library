module lagrange_poly

    implicit none

    public

    contains

    SUBROUTINE lgngr(p,dp,ddp,x,y,nx,ny,type,drctv)
        !***begin prologue     lgngr
        !***date written       940504   (yymmdd)
        !***revision date               (yymmdd)
        !***keywords
        !***author             schneider, barry (nsf)
        !***source             %W% %G%
        !***purpose            lagrange polynomials at arbitrary points.
        !***description
        !***
        !
        !
        !***references
        !
        !***routines called
        !
        !***end prologue       lgngr
        !
        IMPLICIT NONE
        REAL*8, DIMENSION (ny,nx)           :: p
        REAL*8, DIMENSION (ny,nx)           :: dp
        REAL*8, DIMENSION (ny,nx)           :: ddp
        REAL*8, DIMENSION(nx)               :: x
        REAL*8, DIMENSION(ny)               :: y
        REAL*8, DIMENSION(:), ALLOCATABLE   :: xt
        REAL*8, DIMENSION(:), ALLOCATABLE   :: yt
        REAL*8                              :: sn
        REAL*8                              :: ssn
        REAL*8                              :: fac
        CHARACTER (LEN = *), OPTIONAL       :: drctv
        CHARACTER (LEN = *), OPTIONAL       :: type
        INTEGER                             :: nx
        INTEGER                             :: ny
        INTEGER                             :: i
        INTEGER                             :: j
        INTEGER                             :: k
        INTEGER                             :: first
        INTEGER                             :: second
        INTEGER                             :: zerfac
        
        !
        !     generate polynomials and derivatives with respect to x
        !
        p(:,:) = 1.d0
        IF (present(type) ) THEN
           ALLOCATE(xt(nx),yt(ny))
           xt(:) = x(:)
           yt(:) = y(:)
           x(:) = x(:) * x(:)
           y(:) = y(:) * y(:)
        END IF
        DO i=1,ny
          zerfac = 0
          DO j=1,nx
             fac =  y(i) - x(j)
             IF(abs(fac) <= 1.d-10) THEN
                zerfac = j
             ENDIF
          END DO
          DO j=1,nx
             DO k = 1, j-1
                p(i,j) = p(i,j) * ( y(i) - x(k) )   &
                                / ( x(j) - x(k) )
             END DO
             DO k=j+1,nx
                p(i,j) = p(i,j) * ( y(i) - x(k) )   &
                                / ( x(j) - x(k) )
             END DO
             IF(present(drctv) ) THEN
                 IF ( abs(p(i,j)) > 1.d-10) THEN
                      sn = 0.d0
                      ssn = 0.d0
                      DO k=1,j-1
                        !if(i /= k) then
                         fac = 1.d0/( y(i) - x(k) )
                         sn = sn + fac
                         ssn = ssn + fac*fac
                        !endif
                      END DO
                      DO k=j+1,nx
                        !if(i /= k) then
                         fac = 1.d0/( y(i) - x(k) )
                         sn = sn + fac
                         ssn = ssn + fac*fac
                        !end if
                      END DO
                      dp(i,j) = sn*p(i,j)
                      ddp(i,j) = sn*dp(i,j) - ssn*p(i,j)
                 ELSE
                      first=j
                      second=zerfac
                      IF(first > second) THEN
                         first=zerfac
                         second=j
                      END IF
                      sn = 1.d0
                      ssn = 0.d0
                      DO k=1,first-1
                         fac = 1.d0/( x(j) - x(k) )
                         sn = sn*fac*( y(i) - x(k) )
                         ssn = ssn + 1.d0/(y(i) - x(k))
                      END DO
                      DO k=first+1,second-1
                         fac = 1.d0/( x(j) - x(k) )
                         sn = sn*fac*( y(i) - x(k) )
                         ssn = ssn + 1.d0/( y(i) - x(k) )
                      END DO
                      DO k=second+1,nx
                         fac = 1.d0/( x(j) - x(k) )
                         sn = sn*fac*( y(i) - x(k) )
                         ssn = ssn + 1.d0/( y(i) - x(k) )
                      END DO
                      dp(i,j) = sn/( x(j) - x(zerfac) )
                      ddp(i,j) = 2.d0*ssn*dp(i,j)
                 END IF
             END IF
        !
          END DO
        END DO
        !
        IF (present(type)) THEN
           DO i=1,ny
              ddp(i,:) = 2.d0*dp(i,:) + 4.d0 * yt(i) * yt(i) * ddp(i,:)
              dp(i,:) = 2.d0 * yt(i) * dp(i,:)
           END DO
           x(:) = xt(:)
           y(:) = yt(:)
           DEALLOCATE(xt,yt)
        !
        END IF
        
        END SUBROUTINE Lgngr
        
        SUBROUTINE cpoly(cp,dcp,ddcp,pt,n)
        IMPLICIT NONE
        INTEGER                             :: n
        REAL(8), DIMENSION(:,:)           :: cp
        REAL(8), DIMENSION(:,:)           :: dcp
        REAL(8), DIMENSION(:,:)           :: ddcp
        REAL(8), DIMENSION(:)             :: pt
        
        call lgngr(cp,dcp,ddcp,pt,pt,n,n,drctv='on')
        
        end subroutine cpoly

end module