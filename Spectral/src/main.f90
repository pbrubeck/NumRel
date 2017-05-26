program main
    implicit none
    integer :: m,n
    parameter (m=64, n=64)

    !call poisson(n,[1.0d0,1.0d0],[1.0d0,-1.0d0],[1.0d0,-1.0d0])

    call elliptic(m,n)

    print *, 'Program terminated.'
end program main
