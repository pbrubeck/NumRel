program main
    implicit none
    !call poisson(64,[1.0d0,1.0d0],[1.0d0,-1.0d0],[1.0d0,-1.0d0]);
    call elliptic(64,64);
    print *, 'Program terminated.'
end program main
