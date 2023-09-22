!---------------------------------------------------------
! NIST DIV 771
!---------------------------------------------------------
!
! MODULE: matmult
!
! DESCRIPTION:
!> Various generic matrix multiplication routines
!
!---------------------------------------------------------
module matmult
  implicit none

  public mat_mat_mul, mat_vec_mul, tmat_vec_mul, tmat_mat_mul
  private
  real(4), parameter :: s_one  = 1.0
  real(4), parameter :: s_zero = 0.0

  real(8), parameter :: d_one  = 1.d0
  real(8), parameter :: d_zero = 0.d0

  complex(4), parameter :: c_one  = 1.0
  complex(4), parameter :: c_zero = 0.0

  complex(8), parameter :: z_one  = 1.d0
  complex(8), parameter :: z_zero = 0.d0



  !--------------------------------------------------------
  !> @author
  !> Henry Schmale
  !
  ! DESCRIPTION:
  !> generic matrix matrix multiplication interface
  !--------------------------------------------------------
  interface mat_mat_mul
    procedure sfull_mat_mat_mul
    procedure sfull_mat_sym_drealtridiag_mat_mul
    procedure dfull_mat_mat_mul
    procedure dfull_mat_sym_drealtridiag_mat_mul
    procedure cfull_mat_mat_mul
    procedure cfull_mat_sym_drealtridiag_mat_mul
    procedure zfull_mat_mat_mul
    procedure zfull_mat_sym_drealtridiag_mat_mul
  end interface mat_mat_mul

  !--------------------------------------------------------
  !> @author
  !> Henry Schmale
  !
  ! DESCRIPTION:
  !> generic transpose matrix matrix multiplication interface
  !--------------------------------------------------------
  interface tmat_mat_mul
    procedure sfull_tmat_mat_mul
    procedure dfull_tmat_mat_mul
    procedure cfull_tmat_mat_mul
    procedure zfull_tmat_mat_mul
  end interface tmat_mat_mul


  !--------------------------------------------------------
  !> @author
  !> Henry Schmale
  !
  ! DESCRIPTION:
  !> generic matrix vector multiplication interface.
  !> Provides access for many different types
  !--------------------------------------------------------
  interface mat_vec_mul
    procedure sfull_mat_vec_mul
    procedure dfull_mat_vec_mul
    procedure cfull_mat_vec_mul
    procedure zfull_mat_vec_mul
    procedure dfull_mat_zvec_mul
    procedure sym_drealtridiag_mat_dreal_vec_multp
    procedure sym_drealtridiag_mat_dcomplex_vec_multp
    procedure zfull_mat_zvec_mul
    procedure zfull_mat_dvec_mul
    procedure ztridiagonal_mat_zvec_mul
    procedure zmat_band_sym_zvec_mul
    ! procedure ctridiagonal_mat_cvec_mul
    ! procedure ztridiagonal_mat_dvec_mul
    ! procedure ctridiagonal_mat_svec_mul
  end interface mat_vec_mul


  !--------------------------------------------------------
  !> @author
  !> Henry Schmale
  !
  ! DESCRIPTION:
  !> generic transpose matrix vector multiplication interface.
  !> Provides access for many different types
  !--------------------------------------------------------
  interface tmat_vec_mul
    procedure sfull_tmat_vec_mul
    procedure dfull_tmat_vec_mul
    procedure cfull_tmat_vec_mul
    procedure zfull_tmat_vec_mul
    procedure dfull_tmat_zvec_mul
  end interface tmat_vec_mul

  !********************************************
  !********************************************
  !**************************************************************
  !generic, public :: mat_mat_mul  => sf, df,cf, zf
  !generic, public :: tmat_mat_mul => stf, dtf,ctf, ztf
  !generic, public :: mat_vec_mul  => smv, dmv,cmv, zmv, dmzv,dtrimdv, &
  !  dtrimzv, zmdv, zmzv
  !generic, public :: tmat_vec_mul => stmv, dtmv,ctmv, ztmv,dtmzv
  !generic, public :: sym_tridi_mat_vec_mul => dtrimdv,dtrimzv

contains
  !********************************************************
  subroutine sfull_mat_mat_mul( a, b, c)
    real(kind=4), dimension(:,:)                :: a
    real(kind=4), dimension(:,:)                :: b
    real(kind=4), dimension(:,:)                :: c

    integer                                   :: m
    integer                                   :: n
    integer                                   :: k
    integer                                   :: lda
    integer                                   :: ldb
    integer                                   :: ldc

    lda = size(a,1)
    ldb = size(b,1)
    ldc = size(c,1)
    !rows of a and c
    m   = lda
    if(m /= ldc) stop 'wrong matrix dimensions'
    !columns of b and columns of c
    n   = size(b,2)
    if(n /= size(c,2)) stop 'wrong matrix dimensions'
    !rows of b columns of a
    k   = ldb
    if(k /= size(a,2)) stop 'wrong matrix dimensions'

    print*, 'single precision'

    call sgemm ('n', 'n', m, n, k, s_one, a, lda, b, ldb, s_zero, c, ldc)

  end subroutine sfull_mat_mat_mul
  !*********************************************************************
  subroutine dfull_mat_mat_mul( a, b, c)
    real(kind=8), dimension(:,:)                :: a
    real(kind=8), dimension(:,:)                :: b
    real(kind=8), dimension(:,:)                :: c

    integer                                     :: m
    integer                                   :: n
    integer                                   :: k
    integer                                   :: lda
    integer                                   :: ldb
    integer                                   :: ldc

    lda = size(a,1)
    ldb = size(b,1)
    ldc = size(c,1)
    !rows of a and c
    m   = lda
    if(m /= ldc) stop 'wrong matrix dimensions'
    !columns of b and columns of c
    n   = size(b,2)
    if(n /= size(c,2)) stop 'wrong matrix dimensions'
    !rows of b columns of a
    k   = ldb
    if(k /= size(a,2)) stop 'wrong matrix dimensions'

    !print*, 'double precision'

    call dgemm ('n', 'n', m, n, k, d_one, a, lda, b, ldb, d_zero, c, ldc)

  end subroutine dfull_mat_mat_mul

  subroutine cfull_mat_mat_mul( a, b, c)
    complex(kind=4), dimension(:,:)                :: a
    complex(kind=4), dimension(:,:)                :: b
    complex(kind=4), dimension(:,:)                :: c
    integer                                     :: m
    integer                                   :: n
    integer                                   :: k
    integer                                   :: lda
    integer                                   :: ldb
    integer                                   :: ldc

    lda = size(a,1)
    ldb = size(b,1)
    ldc = size(c,1)
    !rows of a and c
    m   = lda
    if(m /= ldc) stop 'wrong matrix dimensions'
    !columns of b and columns of c
    n   = size(b,2)
    if(n /= size(c,2)) stop 'wrong matrix dimensions'
    !rows of b columns of a
    k   = ldb
    if(k /= size(a,2)) stop 'wrong matrix dimensions'

    !print*, 'complex single precision'

    call cgemm ('n', 'n', m, n, k, c_one, a, lda, b, ldb, c_zero, c, ldc)

  end subroutine cfull_mat_mat_mul

  subroutine zfull_mat_mat_mul( a, b, c)
    complex(kind=8), dimension(:,:)                :: a
    complex(kind=8), dimension(:,:)                :: b
    complex(kind=8), dimension(:,:)                :: c
    integer                                     :: m
    integer                                   :: n
    integer                                   :: k
    integer                                   :: lda
    integer                                   :: ldb
    integer                                   :: ldc

    lda = size(a,1)
    ldb = size(b,1)
    ldc = size(c,1)
    !rows of a and c
    m   = lda
    if(m /= ldc) stop 'wrong matrix dimensions'
    !columns of b and columns of c
    n   = size(b,2)
    if(n /= size(c,2)) stop 'wrong matrix dimensions'
    !rows of b columns of a
    k   = ldb
    if(k /= size(a,2)) stop 'wrong matrix dimensions'

    !print*, 'complex double precision'

    call zgemm ('n', 'n', m, n, k, z_one, a, lda, b, ldb, z_zero, c, ldc)

  end subroutine zfull_mat_mat_mul
  !************************************************************************
  !************************************************************************
  subroutine sfull_tmat_mat_mul( a, b, c)
    real(kind=4), dimension(:,:)                :: a
    real(kind=4), dimension(:,:)                :: b
    real(kind=4), dimension(:,:)                :: c

    integer                                     :: m
    integer                                   :: n
    integer                                   :: k
    integer                                   :: lda
    integer                                   :: ldb
    integer                                   :: ldc

    lda = size(a,1)
    ldb = size(b,1)
    ldc = size(c,1)
    !rows of a and c
    m   = lda
    if(m /= ldc) stop 'wrong matrix dimensions'
    !columns of b and columns of c
    n   = size(b,2)
    if(n /= size(c,2)) stop 'wrong matrix dimensions'
    !rows of b columns of a
    k   = ldb
    if(k /= size(a,2)) stop 'wrong matrix dimensions'

   ! print*, 'single precision'

    call sgemm ('t', 'n', m, n, k, s_one, a, lda, b, ldb, s_zero, c, ldc)

  end subroutine sfull_tmat_mat_mul

  subroutine dfull_tmat_mat_mul( a, b, c)
    real(kind=8), dimension(:,:)                :: a
    real(kind=8), dimension(:,:)                :: b
    real(kind=8), dimension(:,:)                :: c

    integer                                     :: m
    integer                                   :: n
    integer                                   :: k
    integer                                   :: lda
    integer                                   :: ldb
    integer                                   :: ldc

    lda = size(a,1)
    ldb = size(b,1)
    ldc = size(c,1)
    !rows of a and c
    m   = lda
    if(m /= ldc) stop 'wrong matrix dimensions'
    !columns of b and columns of c
    n   = size(b,2)
    if(n /= size(c,2)) stop 'wrong matrix dimensions'
    !rows of b columns of a
    k   = ldb
    if(k /= size(a,2)) stop 'wrong matrix dimensions'

    !print*, 'double precision'

    call dgemm ('t', 'n', m, n, k, d_one, a, lda, b, ldb, d_zero, c, ldc)

  end subroutine dfull_tmat_mat_mul

  subroutine cfull_tmat_mat_mul( a, b, c)
    complex(kind=4), dimension(:,:)                :: a
    complex(kind=4), dimension(:,:)                :: b
    complex(kind=4), dimension(:,:)                :: c
    integer                                     :: m
    integer                                   :: n
    integer                                   :: k
    integer                                   :: lda
    integer                                   :: ldb
    integer                                   :: ldc

    lda = size(a,1)
    ldb = size(b,1)
    ldc = size(c,1)
    !rows of a and c
    m   = lda
    if(m /= ldc) stop 'wrong matrix dimensions'
    !columns of b and columns of c
    n   = size(b,2)
    if(n /= size(c,2)) stop 'wrong matrix dimensions'
    !rows of b columns of a
    k   = ldb
    if(k /= size(a,2)) stop 'wrong matrix dimensions'

    !print*, 'complex single precision'

    call cgemm ('t', 'n', m, n, k, c_one, a, lda, b, ldb, c_zero, c, ldc)

  end subroutine cfull_tmat_mat_mul

  subroutine zfull_tmat_mat_mul( a, b, c)
    complex(kind=8), dimension(:,:)                :: a
    complex(kind=8), dimension(:,:)                :: b
    complex(kind=8), dimension(:,:)                :: c
    integer                                     :: m
    integer                                   :: n
    integer                                   :: k
    integer                                   :: lda
    integer                                   :: ldb
    integer                                   :: ldc

    lda = size(a,1)
    ldb = size(b,1)
    ldc = size(c,1)
    !rows of a and c
    m   = lda
    if(m /= ldc) stop 'wrong matrix dimensions'
    !columns of b and columns of c
    n   = size(b,2)
    if(n /= size(c,2)) stop 'wrong matrix dimensions'
    !rows of b columns of a
    k   = ldb
    if(k /= size(a,2)) stop 'wrong matrix dimensions'

    !print*, 'complex double precision'

    call zgemm ('t', 'n', m, n, k, z_one, a, lda, b, ldb, z_zero, c, ldc)

  end subroutine zfull_tmat_mat_mul
  !*********************************************************************
  !*********************************************************************
  !********************************************************
  subroutine sfull_mat_vec_mul( a, b, c, m, n)
    real(kind=4), dimension(:,:)             :: a
    real(kind=4), dimension(:)               :: b
    real(kind=4), dimension(:)               :: c
    integer                                  :: m
    integer                                :: n

    integer                                :: lda

    lda = size(a,1)



    call sgemv ('n', m, n, s_one, a, lda, b, 1, s_zero, c, 1)

  end subroutine sfull_mat_vec_mul
  !*********************************************************************
  subroutine dfull_mat_vec_mul( a, b, c, m, n)
    real(kind=8), dimension(:,:)             :: a
    real(kind=8), dimension(:)               :: b
    real(kind=8), dimension(:)               :: c
    integer                                  :: m
    integer                                :: n

    integer                                :: lda

    lda = size(a,1)



    call dgemv ('n', m, n, d_one, a, lda, b, 1, d_zero, c, 1)

  end subroutine dfull_mat_vec_mul
  !*****************************************************************
  subroutine cfull_mat_vec_mul( a, b, c, m, n)
    complex(kind=4), dimension(:,:)             :: a
    complex(kind=4), dimension(:)               :: b
    complex(kind=4), dimension(:)               :: c
    integer                                     :: m
    integer                                   :: n

    integer                                   :: lda

    lda = size(a,1)



    call cgemv ('n', m, n, c_one, a, lda, b, 1, c_zero, c, 1)

  end subroutine cfull_mat_vec_mul
  !****************************************************************
  subroutine zfull_mat_vec_mul( a, b, c, m, n)
    complex(kind=8), dimension(:,:)             :: a
    complex(kind=8), dimension(:)               :: b
    complex(kind=8), dimension(:)               :: c
    integer                                     :: m
    integer                                   :: n

    integer                                   :: lda

    lda = size(a,1)



    call zgemv ('n', m, n, z_one, a, lda, b, 1, z_zero, c, 1)

  end subroutine zfull_mat_vec_mul
  !********************************************************
  subroutine sfull_tmat_vec_mul( a, b, c, m, n)
    real(kind=4), dimension(:,:)             :: a
    real(kind=4), dimension(:)               :: b
    real(kind=4), dimension(:)               :: c
    integer                                  :: m
    integer                                :: n

    integer                                :: lda

    lda = size(a,1)



    call sgemv ('t', m, n, s_one, a, lda, b, 1, s_zero, c, 1)

  end subroutine sfull_tmat_vec_mul
  !*********************************************************************
  subroutine dfull_tmat_vec_mul( a, b, c, m, n)
    real(kind=8), dimension(:,:)             :: a
    real(kind=8), dimension(:)               :: b
    real(kind=8), dimension(:)               :: c
    integer                                  :: m
    integer                                :: n

    integer                                :: lda

    lda = size(a,1)



    call dgemv ('t', m, n, d_one, a, lda, b, 1, d_zero, c, 1)

  end subroutine dfull_tmat_vec_mul
  !*****************************************************************
  subroutine cfull_tmat_vec_mul( a, b, c, m, n)
    complex(kind=4), dimension(:,:)             :: a
    complex(kind=4), dimension(:)               :: b
    complex(kind=4), dimension(:)               :: c
    integer                                     :: m
    integer                                   :: n

    integer                                   :: lda

    lda = size(a,1)



    call cgemv ('c', m, n, c_one, a, lda, b, 1, c_zero, c, 1)

  end subroutine cfull_tmat_vec_mul
  !****************************************************************
  subroutine zfull_tmat_vec_mul( a, b, c, m, n)
    complex(kind=8), dimension(:,:)             :: a
    complex(kind=8), dimension(:)               :: b
    complex(kind=8), dimension(:)               :: c
    integer                                     :: m
    integer                               :: n
    integer                               :: lda

    lda = size(a,1)

    ! print*,m, n, z_one, a, lda, b, 1, z_zero, c, 1

    call zgemv ('c', m, n, z_one, a, lda, b, 1, z_zero, c, 1)
    !print*,c(1),c(2)
  end subroutine zfull_tmat_vec_mul

  !****************************************************************
  subroutine dfull_mat_zvec_mul( a, b, c, m, n)
    real   (kind=8), dimension(:,:)             :: a
    complex(kind=8), dimension(:)               :: b
    complex(kind=8), dimension(:)               :: c
    integer                                     :: m
    integer                               :: n
    integer                               :: i,j

    c(:) = z_zero

    do j = 1, n
      do i = 1, m

        c(i) = c(i) + a(i,j) * b(j)

      end do
    end do
  end subroutine dfull_mat_zvec_mul

  !*******************************************************************************

  subroutine zfull_mat_dvec_mul( a, b, c)
    complex(kind=8), dimension(:,:)             :: a
    real(kind=8), dimension(:)                  :: b
    complex(kind=8), dimension(:)               :: c
    integer                                     :: m
    integer                               :: n
    integer                               :: i,j

    m = size(a, 1)
    n = size(a, 2)

    c(:) = z_zero

    do j = 1, n
      do i = 1, m
        c(i) = c(i) + a(i,j) * b(j)
      end do
    end do

  end subroutine zfull_mat_dvec_mul

  !*******************************************************************************

  subroutine zfull_mat_zvec_mul( a, b, c)
    complex(kind=8), dimension(:,:)             :: a
    complex(kind=8), dimension(:)               :: b
    complex(kind=8), dimension(:)               :: c
    integer                                     :: m
    integer                                     :: n
    integer                                     :: i,j

    m = size(a, 1)
    n = size(a, 2)

    c(:) = z_zero

    do j = 1, n
      do i = 1, m
        c(i) = c(i) + a(i,j) * b(j)
      end do
    end do

  end subroutine zfull_mat_zvec_mul

  !******************************************************************************************************************


  subroutine dfull_tmat_zvec_mul( a, b, c, m, n)
    real   (kind=8), dimension(:,:)             :: a
    complex(kind=8), dimension(:)               :: b
    complex(kind=8), dimension(:)               :: c
    integer                                     :: m
    integer                               :: n
    integer                               :: i,j

    c(:) = z_zero

    do j = 1, n
      do i = 1, m

        c(i) = c(i) + a(j,i) * b(j)

      end do
    end do



  end subroutine dfull_tmat_zvec_mul

  !******************************************************************************************************************
  ! matrix (tridiagonal and symmetric) vector multiply
  ! takes: real(8) : d =>diagonal of the matrix (ndim)
  !        real(8) : e =>offdiagonal of the matrix (ndim-1)
  !        real(8) : v =>vector (ndim x 1)
  !        integer : ndim => size of the (diagonal of the) matrix
  ! returns: real(8) : ans => ndim x 1 vector
  !*******************************************************************************************************************
  subroutine sym_drealtridiag_mat_dreal_vec_multp(tdmat,v,ans)
    use eigensystem, only : tridiagonal_mat
    type(tridiagonal_mat),intent(in) :: tdmat
    integer                          :: i, ndim
    real(8),intent(in)               :: v(:)
    real(8),intent(out)              :: ans(:)

    ndim = size(v)

    ans(1) =  tdmat%diagonal(1) * v(1)  +  tdmat%offdiagonal (1) * v(2)
    do i = 2, ndim-1
      ans(i) =  tdmat%diagonal(i) * v(i) + tdmat%offdiagonal(i-1) * v(i-1) +  tdmat%offdiagonal (i) * v(i+1)
    end do
    ans(ndim) =  tdmat%diagonal(ndim) * v(ndim) + tdmat%offdiagonal(ndim-1) * v(ndim-1)

  end subroutine sym_drealtridiag_mat_dreal_vec_multp
  !*******************************************************************************************************************
  ! matrix (tridiagonal and symmetric) vector multiply
  ! takes: real(8) : d =>diagonal of the matrix (ndim)
  !        real(8) : e =>offdiagonal of the matrix (ndim-1)
  !     complex(8) : v =>vector (ndim x 1)
  !        integer : ndim => size of the (diagonal of the) matrix
  ! returns: complex(8) : ans => ndim x 1 vector
  !*******************************************************************************************************************
  subroutine sym_drealtridiag_mat_dcomplex_vec_multp(tdmat,v,ans)
    use eigensystem, only : tridiagonal_mat
    type(tridiagonal_mat),intent(in) :: tdmat
    integer                          :: i, ndim
    complex(8),intent(in)            :: v(:)
    complex(8),intent(out)           :: ans(:)

    ndim = size(v)

    ans(1) =  tdmat%diagonal(1) * v(1)  +  tdmat%offdiagonal (1) * v(2)
    do i = 2, ndim-1
      ans(i) =  tdmat%diagonal(i) * v(i) + tdmat%offdiagonal(i-1) * v(i-1) +  tdmat%offdiagonal (i) * v(i+1)
    end do
    ans(ndim) =  tdmat%diagonal(ndim) * v(ndim) + tdmat%offdiagonal(ndim-1) * v(ndim-1)

  end subroutine sym_drealtridiag_mat_dcomplex_vec_multp
  !************************************************************************************************************************

  subroutine sfull_mat_sym_drealtridiag_mat_mul(mat, sym_tdmat,ans)
  use eigensystem, only : tridiagonal_mat
    type(tridiagonal_mat),intent(in) :: sym_tdmat
    integer                          :: i,j, ndim1,ndim2
    real(4),intent(in)            :: mat(:,:)
    real(4),intent(out)           :: ans(:,:)

    ndim1 = size(mat,1)
    ndim2 = size(mat,2)


    do i = 1, ndim1
       ans(i,1) = sym_tdmat%diagonal(1)*mat(1,1) + sym_tdmat%offdiagonal(1)*mat(i,2)
       do j = 2, ndim2-1
          ans(i,j) = sym_tdmat%diagonal(j)*mat(i,j) + sym_tdmat%offdiagonal(j-1)*mat(i,j-1) + &
                       sym_tdmat%offdiagonal(j)*mat(i,j+1)
       end do
       ans(i,ndim2) = sym_tdmat%diagonal(ndim2)*mat(i,ndim2) + sym_tdmat%offdiagonal(ndim2-1)*mat(i,ndim2-1)
    end do

  end subroutine sfull_mat_sym_drealtridiag_mat_mul
  !************************************************************************************************************************

   !************************************************************************************************************************

  subroutine dfull_mat_sym_drealtridiag_mat_mul(mat, sym_tdmat,ans)
  use eigensystem, only : tridiagonal_mat
    type(tridiagonal_mat),intent(in) :: sym_tdmat
    integer                          :: i,j, ndim1,ndim2
    real(8),intent(in)            :: mat(:,:)
    real(8),intent(out)           :: ans(:,:)

    ndim1 = size(mat,1)
    ndim2 = size(mat,2)


    do i = 1, ndim1
       ans(i,1) = sym_tdmat%diagonal(1)*mat(1,1) + sym_tdmat%offdiagonal(1)*mat(i,2)
       do j = 2, ndim2-1
          ans(i,j) = sym_tdmat%diagonal(j)*mat(i,j) + sym_tdmat%offdiagonal(j-1)*mat(i,j-1) + &
                       sym_tdmat%offdiagonal(j)*mat(i,j+1)
       end do
       ans(i,ndim2) = sym_tdmat%diagonal(ndim2)*mat(i,ndim2) + sym_tdmat%offdiagonal(ndim2-1)*mat(i,ndim2-1)
    end do

  end subroutine dfull_mat_sym_drealtridiag_mat_mul
  !************************************************************************************************************************

  !************************************************************************************************************************

  subroutine cfull_mat_sym_drealtridiag_mat_mul(mat, sym_tdmat,ans)
  use eigensystem, only : tridiagonal_mat
    type(tridiagonal_mat),intent(in) :: sym_tdmat
    integer                          :: i,j, ndim1,ndim2
    complex(4),intent(in)            :: mat(:,:)
    complex(4),intent(out)           :: ans(:,:)

    ndim1 = size(mat,1)
    ndim2 = size(mat,2)


    do i = 1, ndim1
       ans(i,1) = sym_tdmat%diagonal(1)*mat(1,1) + sym_tdmat%offdiagonal(1)*mat(i,2)
       do j = 2, ndim2-1
          ans(i,j) = sym_tdmat%diagonal(j)*mat(i,j) + sym_tdmat%offdiagonal(j-1)*mat(i,j-1) + &
                       sym_tdmat%offdiagonal(j)*mat(i,j+1)
       end do
       ans(i,ndim2) = sym_tdmat%diagonal(ndim2)*mat(i,ndim2) + sym_tdmat%offdiagonal(ndim2-1)*mat(i,ndim2-1)
    end do

  end subroutine cfull_mat_sym_drealtridiag_mat_mul
  !************************************************************************************************************************

  !************************************************************************************************************************

  subroutine zfull_mat_sym_drealtridiag_mat_mul(mat, sym_tdmat,ans)
  use eigensystem, only : tridiagonal_mat
    type(tridiagonal_mat),intent(in) :: sym_tdmat
    integer                          :: i,j, ndim1,ndim2
    complex(8),intent(in)            :: mat(:,:)
    complex(8),intent(out)           :: ans(:,:)

    ndim1 = size(mat,1)
    ndim2 = size(mat,2)


    do i = 1, ndim1
       ans(i,1) = sym_tdmat%diagonal(1)*mat(1,1) + sym_tdmat%offdiagonal(1)*mat(i,2)
       do j = 2, ndim2-1
          ans(i,j) = sym_tdmat%diagonal(j)*mat(i,j) + sym_tdmat%offdiagonal(j-1)*mat(i,j-1) + &
                       sym_tdmat%offdiagonal(j)*mat(i,j+1)
       end do
       ans(i,ndim2) = sym_tdmat%diagonal(ndim2)*mat(i,ndim2) + sym_tdmat%offdiagonal(ndim2-1)*mat(i,ndim2-1)
    end do

  end subroutine zfull_mat_sym_drealtridiag_mat_mul
  !************************************************************************************************************************
  subroutine ztridiagonal_mat_zvec_mul(lmat,dmat,umat,vec,ans)

    complex(8)             ,intent(in) :: dmat(:), lmat(:), umat(:)
    complex(8)             ,intent(in) :: vec(:)
    complex(8)             ,intent(out):: ans(:)
    integer                            :: i, ndim

    ndim = size(dmat)

    ans(1) = dmat(1)*vec(1) + umat(1) * vec(2)
    do i = 2, ndim - 1

      ans(i) = dmat(i)*vec(i)
      ans(i) = ans(i) + umat(i) * vec(i+1)
      ans(i) = ans(i) + lmat(i-1) * vec(i-1)

    end do
    ans(ndim) = dmat(ndim)*vec(ndim) + lmat(ndim-1) * vec(ndim-1)

  end subroutine ztridiagonal_mat_zvec_mul

  !//////////////////////////////////////////////////////////////////////////////

  !!*******************************************************************
  !! symetric complex matrix vector multiplier
  !!********************************************************************
  subroutine zmat_band_sym_zvec_mul(zmatd,zmato,vec, ans)
    integer    :: i, j
    integer    :: b, m_size
    complex(8) :: zmatd(:), zmato(:,:)
    complex(8) :: vec(:)
    complex(8) :: ans(:)
    

    b      = size(zmato, 2)
    m_size = size(zmatd)

    ans(1:m_size) = zmatd(1:m_size) * vec(1:m_size)

    do i = 1, b
      ans(1:m_size) = ans(1:m_size) + &
        [zmato(1:m_size-i,i) * vec(i+1:m_size),(z_zero, j=1,i)]
      ans(1:m_size) = ans(1:m_size) + &
        [(z_zero, j=1,i),zmato(1:m_size-i,i) * vec(1:m_size-i)]
    end do

  end subroutine zmat_band_sym_zvec_mul
  !//////////////////////////////////////////////////////////////////////////

end module matmult
