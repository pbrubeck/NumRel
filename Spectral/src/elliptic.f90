subroutine elliptic(m,n)
    implicit none
    integer, intent(in) ::  m, n
    double precision :: V1(m,m-2), L1(m), V2(n,n-2), L2(n)
    double precision :: alph(2,2), beta(2,2)

    alph=reshape([1,1,0,0],[2,2])
    beta=reshape([0,0,1,1],[2,2])

    call helmholtz(m, alph(1,:), beta(1,:), V1, L1)
    call helmholtz(n, alph(2,:), beta(2,:), V2, L2)

    call disp(V1,m,m-2)
    call disp(L1,m,1)
end
