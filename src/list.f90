        module list
        implicit none
        
        type cell
           real(8) :: value
           integer :: indx
           type(cell),pointer :: nextcell
        end type cell
        
        type cell_array
           real(8), dimension(:), allocatable :: array
           integer                            :: indx
           type(cell_array),pointer           :: nextcell
        end type cell_array
        
        
        contains
!****************************************        
        subroutine init_fill_list(first,aElem,data)
        real(8),intent(in)  :: data
        integer             :: AllocStat
        type(cell), pointer :: first, aElem
        
        !Writing the list
           allocate( first, stat = AllocStat )
           first%value = data        
           first%indx  = 1
           aElem => first
        
        end subroutine
!****************************************        
        subroutine fill_list(aElem,data)
        real(8),intent(in)  :: data
        type(cell), pointer :: aElem
                
           allocate(aElem%nextcell)
           aElem%nextcell%value = data
           aElem%nextcell%indx = aElem%indx + 1
           aElem => aElem%nextcell
           
        end subroutine
!****************************************        
        subroutine end_list(aElem)
        type(cell), pointer ::  aElem
        
           nullify(aElem%nextcell)
        
        end subroutine
!****************************************        
        !Now reading back
        
        subroutine read_list(first, aElem, data, length)
        real(8),allocatable :: data(:)
        integer,intent(out) :: length
        integer             :: i
        type(cell), pointer :: first, aElem
          length = aElem%indx
          allocate(data(length))
          aElem => first
          i = 0
          do while(associated(aElem))
            i = i + 1
            data(i) = aElem%value
            aElem => aElem%nextcell
          end do
        end subroutine  
!****************************************
! list of a one dimensional array
!****************************************
!****************************************        
        subroutine init_fill_array_list(first,aElem,data)
        real(8),intent(in),dimension(:)  :: data
        integer                          :: AllocStat
        type(cell_array), pointer        :: first, aElem
                
        !Writing the list
           allocate( first, stat = AllocStat )
           allocate( first%array(size(data)) )
           
           first%array = data        
           first%indx  = 1
           aElem => first
        
        end subroutine
!****************************************        
        subroutine fill_array_list(aElem,data)
        real(8),intent(in),dimension(:)  :: data
        type(cell_array), pointer        :: aElem
                
           allocate(aElem%nextcell)
           !allocate(aElem%array(size(data)))
           !allocate(aElem%nextcell%array(size(data)))
           aElem%nextcell%array = data
           aElem%nextcell%indx = aElem%indx + 1
           aElem => aElem%nextcell
           
        end subroutine
!****************************************        
        subroutine end_array_list(aElem)
        type(cell_array), pointer ::  aElem
        
           nullify(aElem%nextcell)
        
        end subroutine
! !****************************************        
        !Now reading back
        
        subroutine read_array_list(first, aElem, data, length)
        real(8),allocatable       :: data(:,:)
        integer,intent(out)       :: length
        integer                   :: i
        type(cell_array), pointer :: first, aElem
          length = aElem%indx
          allocate(data(size(first%array),length))
          aElem => first
          i = 0
          do while(associated(aElem))
            i = i + 1
            data(:,i) = aElem%array(:)
            aElem => aElem%nextcell
          end do
        end subroutine  
! !****************************************        

        end module
