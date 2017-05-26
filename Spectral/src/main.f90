program main
    implicit none

    !call poisson(n,[1.0d0,1.0d0],[1.0d0,-1.0d0],[1.0d0,-1.0d0])

    call elliptic(32,32)


    print *, 'Program terminated.'
end program main
