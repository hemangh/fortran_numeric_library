module eigensystem

  implicit none
  private
  public tridiagonal_mat, tri_diagonalizer, all_vec, some_vec, &
         nstridiagonal_mat

  !************************************************************************
  interface
    !! Lapack function that determines double precision machine parameters.
    function dlamch(ch)
      character(1)::ch
      real(8) :: DLAMCH
    end function dlamch
  end interface
  !*************************************************************************
  ! this refers to a symmetric tridiagonal matrix
  !*********************************************************************
  type tridiagonal_mat
    integer                       :: mat_size
    real(8), allocatable          :: diagonal   (:)
    real(8), allocatable          :: offdiagonal(:)
  end type
  !*********************************************************************
  ! this refers to a non-symmetric tridiagonal matrix
  !*************************************************************************
  type nstridiagonal_mat
    integer                       :: mat_size
    real(8), allocatable          :: diagonal   (:)
    real(8), allocatable          :: uoffdiagonal(:)
    real(8), allocatable          :: loffdiagonal(:)
  end type
  !*********************************************************************

  !! General tri_diagonal matrix eigensystem solver's parameters
  !! that are shared wheather all or some eigenvectors need to be found.
  type, extends(tridiagonal_mat) :: tri_diagonalizer

    integer :: info
    integer :: ldz
    !allocatable variables
    real(8), allocatable :: WORK( : ), z(:,:), eigenval(:),eigenvec(:,:)

  contains

    procedure :: mat_initialize
    procedure :: allocation
    procedure :: deallocation
    procedure :: eigensystem_solve

  end type
  !****************************************************************************
  !****************************************************************************
  !! Extended tri_diagonal matrix eigensysem solver's parameters for the case
  !! where all eigenvectors are desired
  type, extends (tri_diagonalizer) :: all_vec
    !     .. Scalar Arguments ..
    CHARACTER  ::       compz

  end type
  !*****************************************************************************
  !*****************************************************************************
  !! Extended tri_diagonal matrix eigensystem solver's parameters for the case
  !! where only some eigenvectors are desired.
  type, extends (tri_diagonalizer) :: some_vec
    !     .. Array Arguments ..
    INTEGER,allocatable          ::  IBLOCK( : ), IFAIL( : ), ISPLIT( : ),&
      IWORK( : )
    DOUBLE PRECISION,allocatable ::  W( : )
    !     .. Scalar Arguments ..
    CHARACTER                    ::  ORDER, RANGE
    INTEGER                      ::  M
    INTEGER                      ::  IL, IU, NSPLIT
    DOUBLE PRECISION             ::  ABSTOL, VL, VU, DLAMCH

  end type
  !******************************************************************************

contains

  ! !*********************************************************************
  subroutine mat_initialize(mat,d,e,ldz,dir1,mfind,dir2,il,iu)

    class (tri_diagonalizer) :: mat
    real(8) :: d(:), e(:)
    integer :: n ! is sizeof d
    integer, optional :: ldz, mfind,il,iu
    character(len=1),optional:: dir1,dir2
    ! real*8                   :: dlamch

    n = size(d)

    allocate(mat%diagonal(n),mat%offdiagonal(n-1))
    mat%mat_size = n
    mat%diagonal = d
    mat%offdiagonal = e
    if (present(ldz)) then
      mat%ldz=ldz
    else
      mat%ldz= n
    end if

    select type (mat)
    type is (tri_diagonalizer)

      class is (all_vec)
      ! from lapack dsteqr.f
      ! *  COMPZ   (input) CHARACTER*1
      ! *          = 'N':  Compute eigenvalues only.
      ! *          = 'V':  Compute eigenvalues and eigenvectors of the original
      ! *                  symmetric matrix.  On entry, Z must contain the
      ! *                  orthogonal matrix used to reduce the original matrix
      ! *                  to tridiagonal form.
      ! *          = 'I':  Compute eigenvalues and eigenvectors of the
      ! *                  tridiagonal matrix.  Z is initialized to the identity
      ! *                  matrix.
      if (present(dir1)) then
        mat%compz=dir1
      else
        mat%compz= "I"
      end if
      class is (some_vec)
      mat%abstol=2*DLAMCH('S')
      if (present(dir1)) then
        mat%range=dir1
        if(dir1 == 'I')then
          if( present(il) .and. present(iu))then
            mat%il = il
            mat%iu = iu
          else
            print*,'Forgot to initialize il and iu for some_vec&
              eigenvals'
            stop
          end if
        end if
      else
        mat%range= "A"
      end if
      if (present(mfind)) then
        mat%m =mfind
      else
        mat%m = n
      end if
      if (present(dir2)) then
        mat%order=dir2
      else
        mat%order= "E"
      end if

      class default
      ! ! give error for unexpected/unsupported type
      print*,'allocation: unexpected type for mat object!'
      stop

    end select

    if(size(e) /= n-1)then
      print*,'error with offdiagonal size of array'
      stop
    end if

    call allocation(mat)

  end subroutine mat_initialize
  ! !***********************************************************************
  subroutine allocation(mat)
    class (tri_diagonalizer) :: mat

    select type (mat)
    type is (tri_diagonalizer)

      class is (all_vec)
      allocate(mat%z(mat%ldz,mat%mat_size))
      allocate(mat%work(2*(mat%mat_size-1)))
      allocate(mat%eigenvec(mat%ldz,mat%mat_size))
      allocate(mat%eigenval(mat%mat_size))
      class is (some_vec)
      allocate(mat%z(mat%ldz,mat%m),mat%w(mat%mat_size))
      allocate(mat%work(5*mat%mat_size),&
        mat%ISPLIT(mat%mat_size),mat%IFAIL(mat%mat_size),&
        mat%IBLOCK(mat%mat_size),mat%IWORK(3*mat%mat_size))
      allocate(mat%eigenvec(mat%ldz,mat%m))
      allocate(mat%eigenval(mat%mat_size))

      class default
      ! ! give error for unexpected/unsupported type
      print*,'allocation: unexpected type for mat object!'
      stop
    end select



  end subroutine allocation
  ! !**************************************************************************
  subroutine deallocation(mat)
    class (tri_diagonalizer) :: mat

    select type (mat)
    type is (tri_diagonalizer)

      class is (all_vec)
      deallocate(mat%diagonal,mat%offdiagonal)
      deallocate(mat%z)
      deallocate(mat%work)
      deallocate(mat%eigenvec)
      deallocate(mat%eigenval)
      class is (some_vec)
      deallocate(mat%diagonal,mat%offdiagonal)
      deallocate(mat%z,mat%w)
      deallocate(mat%work,&
        mat%ISPLIT,mat%IFAIL,&
        mat%IBLOCK,mat%IWORK)
      deallocate(mat%eigenvec)
      deallocate(mat%eigenval)

      class default
      ! ! give error for unexpected/unsupported type
      print*,'allocation: unexpected type for mat object!'
      stop
    end select



  end subroutine deallocation

  ! !**************************************************************************
  subroutine eigensystem_solve(mat)

    integer :: m
    class (tri_diagonalizer) :: mat

    select type (mat)
    type is (tri_diagonalizer)
      ! ! no further allocation is needed
      class is (all_vec)
      ! print*, 'inside diagonal:', mat%diagonal

      !print*,'All eigenvalues and eigenvectors are computed'
      call dsteqr (mat%COMPZ, mat%mat_size, mat%diagonal, mat%offdiagonal, &
        mat%Z, mat%LDZ, mat%WORK, mat%INFO)

      if (mat%info /= 0) then
        print*, 'Info =', mat%info
        stop
      else
        mat%eigenval= mat%diagonal
        mat%eigenvec= mat%z
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      class is (some_vec)
      !print*,'All eigenvalues and some eigenvectors are computed'

      call DSTEBZ( mat%RANGE, mat%ORDER, mat%mat_size, mat%VL, mat%VU, &
        mat%IL, mat%IU, mat%ABSTOL, mat%diagonal, mat%offdiagonal,  &
        m, mat%NSPLIT, mat%W, mat%IBLOCK, mat%ISPLIT, mat%WORK, &
        mat%IWORK,mat%INFO )
      if (mat%info /= 0) then
        print*, 'Info =', mat%info
        stop
      else
        mat%eigenval= mat%w
      end if

      call dstein( mat%mat_size, mat%diagonal, mat%offdiagonal, mat%M, mat%W, &
        mat%IBLOCK, mat%ISPLIT, mat%Z, mat%LDZ,mat%WORK, mat%IWORK,&
        mat%IFAIL, mat%INFO )

      if (mat%info /= 0) then
        print*, 'Info =', mat%info
        stop
      else
        mat%eigenvec= mat%z
      end if

      class default
      ! ! give error for unexpected/unsupported type
      print*,'allocation: unexpected type for mat object!'
      stop
    end select



  end subroutine eigensystem_solve
  !*******************************************************************************
end module eigensystem
