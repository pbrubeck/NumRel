subroutine elliptic(m,n)
    implicit none
    integer, intent(in) :: m, n
    integer :: i, j, lwork, info, piv1(m-2), piv2(n-2)
    double precision :: alph(2,2), beta(2,2)
    double precision :: x(m), y(n), bc1(2,n), bc2(m,2), F(m,n)
    double precision :: A1(m,m), D1(m,m), V1(m,m-2), W1(m-2,m-2), L1(m), B1(2,m), N1(m,2)
    double precision :: A2(n,n), D2(n,n), V2(n,n-2), W2(n-2,n-2), L2(n), B2(2,n), N2(n,2)
    double precision :: uu(m,n), u1(m,n), u2(m,n), c1(m,n), c2(m,n), aux1(m,n), aux2(n,m)
    double precision :: Lij, WORK(8*(m+n))
    double precision :: pi
    parameter (pi=4*atan(1.0d0))

    lwork=8*(m+n);

    alph=reshape([1,1,1,1],[2,2]);
    beta=reshape([0,0,0,0],[2,2]);

    call helmholtz(m, x, D1, A1);
    call helmholtz(n, y, D2, A2);
    call setbc(m, alph(1,:), beta(1,:), D1, A1, V1, L1, B1, N1);
    call setbc(n, alph(2,:), beta(2,:), D2, A2, V2, L2, B2, N2);

    bc1(1,:)=0.2*sin(3*pi*y);
    bc1(2,:)=0*y;
    bc2(:,1)=merge(sin(pi*x)**4, 0.0d0, x<0);
    bc2(:,2)=0*x;
    F(:,:)=0.0d0;

    ! C=(B1*B1')*(bc1*B2')+(B1*bc2)*(B2*B2');
    ! uu=N1*sylvester(B1*B1', B2*B2', C)*N2';
    uu(:,:)=0.0d0;

    ! uu=uu+(N1*bc1-uu)*P2'+P1*(bc2*N2'-uu);
    ! P=I-B'/(B*B')*B

    ! u1 = N1*bc1-uu
    ! u2 = bc2*N2'-uu
    ! uu = u+u1*P2+P1*u2

    u1=uu; ! TODO: deep copy
    u2=uu; ! TODO: deep copy
    call dgemm('N', 'N', m, n, 2, 1.0d0, N1, m, bc1, 2, -1.0d0, u1, m);
    call dgemm('N', 'T', m, n, 2, 1.0d0, bc2, m, N2, n, -1.0d0, u2, m);

    ! u1=u1*P2 TODO: complementary projector
    ! u2=P1*u2 TODO: complementary projector
    uu=uu+u1+u2;

    ! Right hand side
    ! F=F-A1*uu-uu*A2';
    call dgemm('N', 'N', m, n, m, -1.0d0, A1, m, uu, m, 1.0d0, F, m);
    call dgemm('N', 'T', m, n, n, -1.0d0, uu, m, A2, n, 1.0d0, F, m);


    ! Inversion via diagonalization
    ! F(2:m-1,2:n-1)=V1(2:m-1,:)\F(2:m-1,2:n-1)/V2(2:n-1,:)';
    W1=V1(2:m-1,:);
    W2=V2(2:n-1,:);
    call dgetrf(m-2, m-2, W1, m-2, piv1, info);
    call dgetrf(n-2, n-2, W2, n-2, piv2, info);
    call dgetrs('N', m-2, m-2, W1, m-2, piv1, F(2:m-1,2:n-1) , m-2, info);
    aux2(2:n-1,2:m-1)=transpose(F(2:m-1,2:n-1));
    call dgetrs('N', n-2, n-2, W2, n-2, piv2, aux2(2:n-1,2:m-1), n-2, info);
    F(2:m-1,2:n-1)=transpose(aux2(2:n-1,2:m-1));
    ! F(2:m-1,2:n-1)=F(2:m-1,2:n-1)./LL;
    do i=2,m-1
        do j=2,n-1
            Lij=L1(i)+L2(j);
            F(i,j)=F(i,j)/Lij;
        enddo
    enddo

    ! uu=uu+V1*F(2:m-1,2:n-1)*V2'
    call dgemm('N', 'N', m, n-2, m-2, 1.0d0, V1, m, F(2:m-1,2:n-1), m-2, 0.0d0, aux1(:,2:n-1), m);
    call dgemm('N', 'T', m, n  , n-2, 1.0d0, aux1(:,2:n-1), m, V2, n, 1.0d0, uu, m);
    call disp(uu,m,n);
end
