subroutine helmholtz(n, x,D,D2)
    implicit none
    integer, intent(in) :: n
    double precision :: x(n), D(n,n), D2(n,n)
    call chebD(n, D, x);
    call dgemm('N', 'N', n, n, n, 1.0d0, D, n, D, n, 0.0d0, D2, n);
end
