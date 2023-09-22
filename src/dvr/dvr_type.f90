module dvr_type
    use,intrinsic :: iso_fortran_env
    use radial_grid_type
    use lagrange_poly
    implicit none

    PUBLIC :: dvr_set, renormalize_dvr_set_array
    PRIVATE

    integer, parameter :: idp = real64
!************************************************************************************************
! dvr_set is an extension of the radial_grid for a contiguous 
!  set in one region
!**************************************************************************************************
type  :: dvr_set
  type(radial_grid), dimension(:), allocatable :: radial_grid
  integer                                      :: num_pts
  real(idp), dimension(:),allocatable          :: norm      ! normalization of lobotto Gauss quadratures in each region dimension(region_num)
  real(idp), dimension(:),allocatable          :: pw         !lagrange polynomial weights (for boundary they are different)
  real(idp), dimension(:,:),allocatable        :: p          !lagrange polynomials of FEDVRS, has dimension (num,num) followed by:
  real(idp), dimension(:,:),allocatable        :: dp         !first derivative
  real(idp), dimension(:,:),allocatable        :: ddp        !second derivative
  logical                                      :: called = .false.
contains
  procedure, private :: construct_dvr_set
  generic, public    :: construct => construct_dvr_set
  procedure, private :: destruct_dvr_set
  generic, public    :: destruct => destruct_dvr_set
end type

contains

! ****************************************************************************************
! input: 
!       this: the dvr class
!       rad_grid: the radial_grid type that has to be already constructed and be legendre
!                 gaussian (e.g. natural grid of the dvrs)
!                 if a different set of points are being used than the natural 
!                 fedvr points (e.g. Gauss-Legendre points) still a dvr set can be formed
!                 however, attention must be paid to the normalization and weights of such
!                 a set (recommended this option not be used)
! ***************************************************************************************************
subroutine construct_dvr_set(this,rad_grid)
    class(dvr_set) :: this
    type(radial_grid), intent(in) :: rad_grid(:)

    integer :: num_pts, num_regions
    integer :: iregion, left_indx, right_indx

    real(idp), allocatable :: p(:,:), dp(:,:), ddp(:,:)

    num_regions = size(rad_grid)

    num_pts = 0

    do iregion = 1, num_regions
     if (rad_grid(iregion)%called .eqv. .false.) stop 'dvr set concstruct called without an initialized grid'
     num_pts = num_pts + rad_grid(iregion)%num_pts
    end do

     num_pts = num_pts - num_regions + 1

     this%num_pts = num_pts
     this%radial_grid = rad_grid
     allocate(this%pw(num_pts), this%norm(num_pts))
     allocate(this%p(num_pts,num_pts), this%dp(num_pts,num_pts), this%ddp(num_pts,num_pts))
     
     !initialize:
     this%p(:,:) = 0._idp
     this%dp(:,:) = 0._idp
     this%ddp(:,:) = 0._idp

     left_indx = 1
     right_indx = rad_grid(1)%num_pts
     this%pw(left_indx:right_indx) = rad_grid(1)%weights
     call cpoly(this%p(left_indx:right_indx, left_indx:right_indx), &
                 this%dp(left_indx:right_indx, left_indx:right_indx), &
                 this%ddp(left_indx:right_indx, left_indx:right_indx), &
                 rad_grid(1)%pts,right_indx)

     do iregion = 2, num_regions
            num_pts = rad_grid(iregion)%num_pts
            left_indx = right_indx
            this%pw(left_indx) = this%pw(left_indx) + rad_grid(iregion)%weights(1)
            right_indx = right_indx + num_pts - 1
            this%pw(left_indx+1:right_indx) = rad_grid(iregion)%weights(2:)
            allocate(p(num_pts,num_pts),dp(num_pts,num_pts),ddp(num_pts,num_pts))
            call cpoly(p,dp,ddp,rad_grid(iregion)%pts,num_pts)
            this%p(left_indx+1:right_indx, left_indx+1:right_indx) = p(2:,2:)
            this%dp(left_indx+1:right_indx, left_indx+1:right_indx) = dp(2:,2:)
            this%dp(left_indx,left_indx) = this%dp(left_indx,left_indx) + dp(1,1)
            this%ddp(left_indx+1:right_indx, left_indx+1:right_indx) = ddp(2:,2:)
            this%ddp(left_indx,left_indx) = this%ddp(left_indx,left_indx) + ddp(1,1)
            deallocate(p,dp,ddp)
     end do    
     this%norm = 1._idp / sqrt(this%pw)
    
    ! !
    this%called=.true.

    
    end subroutine

subroutine destruct_dvr_set(this)
    class (dvr_set) :: this
    ! if (allocated(this%pts) ) then
    !   deallocate(this%pts, this%pw, this%p, this%dp, this%ddp, this%norm)
    ! end if
    ! this%called = .false.
    end subroutine


subroutine renormalize_dvr_set_array (dvr_set_array)   
  
  type(dvr_set) :: dvr_set_array (:)
  ! integer       :: num_regions, i, left_indx

  ! num_regions = size (dvr_set_array) 
  ! if (num_regions > 1) then

 
  !   do i= 1, num_regions-1
       
  !      left_indx = dvr_set_array(i)%num_pts
  !     !  if (dvr_set_array(i)%pts(left_indx) == dvr_set_array(i+1)%pts(1)) then
  !        dvr_set_array(i)%norm(left_indx) = 1._idp / sqrt(dvr_set_array(i)%pw(left_indx) + dvr_set_array(i+1)%pw(1))
  !        dvr_set_array(i+1)%norm(1)       = dvr_set_array(i)%norm(left_indx)
  !        dvr_set_array(i)%pw(left_indx) = dvr_set_array(i)%pw(left_indx) + dvr_set_array(i+1)%pw(1)
  !        dvr_set_array(i+1)%pw(1)       = dvr_set_array(i)%pw(left_indx)
  !     !  end if
 
  !   end do
 
 
  ! end if

 end subroutine 



end module