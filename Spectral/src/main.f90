program main
    implicit none
    !call poisson(64,[1.0d0,1.0d0],[1.0d0,-1.0d0],[1.0d0,-1.0d0]);
    call elliptic(6,6);
    print *, 'Program terminated.'
end program main
