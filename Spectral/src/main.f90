program main
    implicit none
    !call poisson(64,[1.0d0,1.0d0],[1.0d0,-1.0d0],[1.0d0,-1.0d0]);

    call elliptic(128,128);
    print *, 'Program terminated.'
end program main
