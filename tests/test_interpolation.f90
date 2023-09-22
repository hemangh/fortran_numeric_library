program test_interpolation
    use interpolation, only: linear_intp
    implicit none
    
    integer n
    double precision, allocatable :: x(:), y(:)
    double precision ux, uy
    
    n = 2
    allocate(x(n), y(n))
    x = [1.d0,2.d0]
    y = [1.d0,2.d0]
    ux = 1.5d0
    uy = linear_intp (ux, x, y, n)
    if(uy /= 1.5d0) call exit(1)

    ux = 3.d0
    uy = linear_intp (ux, x, y, n)
    

    deallocate(x,y)



end program test_interpolation