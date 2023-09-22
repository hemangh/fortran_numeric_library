program test_list
    use list
    implicit none
    

    type(cell), pointer :: list_init, list_member
    real(8):: harvest = 0.d0
    real(8), allocatable :: data(:)
    integer :: i, length
     call init_fill_list(list_init, list_member, harvest)

    do i = 1, 10000
        call random_number(harvest)
        call fill_list(list_member, harvest)
    end do
    call end_list(list_member)
    length = 10001
    call read_list(list_init, list_member, data, length)
    
end program test_list