module radial_grid_type
    use,intrinsic :: iso_fortran_env
    use gauss_quadrature_mod
    implicit none

    PUBLIC :: radial_grid
    PRIVATE

    integer, parameter :: idp = real64

  !*************************************************************************************************************************************
  !radial type (public) radial 
    ! this type includes components for a radial grid as well as the ability to interpolate
    ! on such a grid
  !*************************************************************************************************************************************
  type radial_grid
    integer                                      :: num_pts
    integer, dimension(2)                        :: fixed_end !determines if the boundaries are included in quadrature (1), default is not (0).
    integer                                      :: kpts
    character(len=32)                            :: ckind     ! type of grid (keywords listed in dvr.f90)
    real(idp), dimension(:),allocatable          :: pts         ! the points of the grid
    real(idp), dimension(:),allocatable          :: weights         ! the weights of the grid points 
    real(idp), dimension(:),allocatable          :: dr3
    real(idp), dimension(2)                      :: endpts
    real(idp)                                    :: alpha
    real(idp)                                    :: beta
    real(idp)                                    :: left_edge
    real(idp)                                    :: right_edge
    logical  , dimension(2)                      :: drop=.false.!determines whether to fix then drop either of the endpoints.
    logical                                      :: called=.false.
  contains
    procedure, private :: construct_radial_type
    procedure, private :: destroy_radial_type
    generic, public :: construct => construct_radial_type
    generic, public :: destruct  =>  destroy_radial_type

end type


contains
! ****************************************************************************************
! input: 
!       this: the radial class
!       num_pts: integer, number of points in the radial region
!       end_pts: real(2), bounderies or edges of the radial region. index 1 -> left boundary
!                index 2 -> right boundary
!       grid_kind: character, describes the kind of grid being constructed
!                  all the different kinds of the grids available are described
!                  in the Library/dvr/dvr.f90 -> function_class
!       fixed_end: integer(2), has values of only 0 and 1, indicates if the grid includes either of
!                  the two end points. index 1 -> left boundary (1: included, 0: not included)
!                  index 2 -> right boundary (1: included 0: not included)
!       fix_and_drop: logical(2) option to drop the first and/or last points in the grid 
!                     it may be usefull, but generally not recommended to drop a point after making 
!                     the grid. 
! ***************************************************************************************************
subroutine construct_radial_type(this,num_pts,end_pts, grid_kind, fixed_end, fix_and_drop)
    class(radial_grid) :: this
    integer, intent(in) :: num_pts
    real(idp), intent(in), dimension(2) :: end_pts
    integer, dimension(2), intent(in)   :: fixed_end
    character (len=*), intent(in) :: grid_kind
    logical, dimension(2), optional :: fix_and_drop
    
    this%num_pts = num_pts
    allocate(this%pts(num_pts), this%weights(num_pts))
    this%left_edge  = end_pts(1)
    this%right_edge = end_pts(2)
    this%fixed_end = fixed_end
    this%ckind = grid_kind

    if(this%fixed_end(1)==1 .and. this%fixed_end(2)==0)then
        !left point is only fixed
        this%kpts=1
        this%endpts(1)= -1.d0
        this%endpts(2)= 1.d0
      elseif(this%fixed_end(1)==0 .and. this%fixed_end(2)==1)then
        this%kpts=1
        this%endpts(1)= 1.d0
        this%endpts(2)= -1.d0
      elseif(this%fixed_end(1)==1 .and. this%fixed_end(2)==1)then
        !both are fixed
        this%kpts=2
        this%endpts(1)=-1.d0
        this%endpts(2)= 1.d0
      else
        !none of the end points are fixed (endpts is irrelevant) or
        this%kpts = 0
      end if

    call gauss_quadrature(trim(this%ckind),this%num_pts, this%kpts,&
                          this%endpts, this%pts, this%weights)
    if (present(fix_and_drop))then
        this%drop = fix_and_drop
    end if    
    this%endpts(1)= this%left_edge
    !
    this%endpts(2)= this%right_edge
    
    ! mapping the points to the two real ends
    call map_to_endpoints(this%pts,this%weights,this%endpts,this%num_pts)
    !
    ! left_edge or right_edge at this point could be different than end_pts
    ! distinction here is that even though the boundaries are set in end_pts
    ! a grid point might be excluded or the first or last grid point might 
    ! not fall on the boundary points. 
    this%left_edge   = this%pts(1)
    this%right_edge  = this%pts(num_pts)

    ! assigning r^2 dr
    this%dr3 = this%weights * this%pts * this%pts
    !
    !
    this%called=.true.

    
    end subroutine

subroutine destroy_radial_type(this)
    class (radial_grid) :: this
    if (allocated(this%pts) ) then
    deallocate(this%pts, this%weights, this%dr3)
    end if
    this%called = .false.
    end subroutine

end module