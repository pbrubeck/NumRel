subroutine elliptic(m,n)
    ! Solves separatable elliptic equations
    use linalg
    use visual
    implicit none
    integer, intent(in) :: m, n
    integer :: i, j, info, piv1(m-2), piv2(n-2)
    double precision :: alph(4), beta(4)
    double precision :: x(m), y(n), bc1(2,n), bc2(m,2), F(m,n)
    double precision :: A1(m,m), D1(m,m), V1(m,m-2), W1(m-2,m-2), L1(m), B1(2,m), N1(m,2)
    double precision :: A2(n,n), D2(n,n), V2(n,n-2), W2(n-2,n-2), L2(n), B2(2,n), N2(n,2)
    double precision :: uu(m,n), u0(2,2), u1(m,n), u2(m,n), BB1(2,2), BB2(2,2)
    double precision :: aux1(m,n), aux2(n,m)
    double precision :: Lij
    double precision, parameter :: pi=4*atan(1.0d0)

    alph=[1.0d0, 0.0d0, 1.0d0, 0.0d0];
    beta=[0.0d0, 1.0d0, 0.0d0, 1.0d0];

    ! Differential operators
    ! TODO: user supplied
    call helmholtz(m, x, D1, A1);
    call helmholtz(n, y, D2, A2);

    ! Find eigenfunctions V, boundary condition operator B and dual N
    ! This is the most computationally expensive part
    call setbc(m, alph(1:2), beta(1:2), D1, A1, V1, L1, B1, N1);
    call setbc(n, alph(3:4), beta(3:4), D2, A2, V2, L2, B2, N2);

    ! Problem data
    ! TODO: user supplied
    ! Boundary values
    bc1(1,:)=0.2*sin(3*pi*y);
    bc1(2,:)=0*y;
    bc2(:,1)=merge(sin(pi*x)**4, 0.0d0, x<0);
    bc2(:,2)=0*x;
    ! Source term
    F(:,:)=0.0d0;

    ! Non-homogenous contribution
    ! Low-rank multiplications, matmul is justified
    ! u0=(B1*B1')*(bc1*B2')+(B1*bc2)*(B2*B2');
    ! uu=N1*sylvester(B1*B1', B2*B2', u0)*N2';
    BB1=matmul(B1, transpose(B1));
    BB2=matmul(B2, transpose(B2));
    u0=matmul(BB1, matmul(bc1, transpose(B2)))+matmul(matmul(B1, bc2), BB2);
    u0=sylvester2(BB1, BB2, u0);
    uu=matmul(matmul(N1,u0), transpose(N2));

    ! uu=uu+(N1*bc1-uu)*P2'+P1*(bc2*N2'-uu);
    ! Complementary projector P=I-B'/(B*B')*B
    ! u1=N1*bc1-uu
    ! u2=bc2*N2'-uu
    ! uu=uu+u1*P2+P1*u2
    u1=uu;
    u2=uu;
    call dgemm('N', 'N', m, n, 2, 1.0d0, N1, m, bc1, 2, -1.0d0, u1, m);
    call dgemm('N', 'T', m, n, 2, 1.0d0, bc2, m, N2, n, -1.0d0, u2, m);
    u1=u1-matmul(matmul(u1, transpose(B2)), matmul(matinv2(BB2), B2));
    u2=u2-matmul(matmul(transpose(B1), matinv2(BB1)), matmul(B1, u2));
    uu=uu+u1+u2;

    W1=V1(2:m-1,:);
    W2=V2(2:n-1,:);
    call dgetrf(m-2, m-2, W1, m-2, piv1, info);
    call dgetrf(n-2, n-2, W2, n-2, piv2, info);

    ! Pre-computation phase finished

    ! Right-hand side
    ! F=F-A1*uu-uu*A2';
    call dgemm('N', 'N', m, n, m, -1.0d0, A1, m, uu, m, 1.0d0, F, m);
    call dgemm('N', 'T', m, n, n, -1.0d0, uu, m, A2, n, 1.0d0, F, m);

    ! Inversion via diagonalization
    ! Analogous to Green's function

    ! First, eigenfunction expansion
    ! F(2:m-1,2:n-1)=V1(2:m-1,:)\F(2:m-1,2:n-1)/V2(2:n-1,:)';
    call dgetrs('N', m-2, n-2, W1, m-2, piv1, F(2:m-1,2:n-1) , m-2, info);
    ! TODO: Fix this transposition
    aux2(2:n-1,2:m-1)=transpose(F(2:m-1,2:n-1));
    call dgetrs('N', n-2, m-2, W2, n-2, piv2, aux2(2:n-1,2:m-1), n-2, info);
    F(2:m-1,2:n-1)=transpose(aux2(2:n-1,2:m-1));

    ! Second, divide by the eigenvalues
    ! F(2:m-1,2:n-1)=F(2:m-1,2:n-1)./LL;
    F(2:m-1,2:n-1)=F(2:m-1,2:n-1)/(spread(L1(2:m-1),2,n-2)+spread(L2(2:n-1),1,m-2));

    ! Third, eigenfunction synthesis
    ! uu=uu+V1*F(2:m-1,2:n-1)*V2'
    call dgemm('N', 'N', m, n-2, m-2, 1.0d0, V1, m, F(2:m-1,2:n-1), m-2, 0.0d0, aux1(:,2:n-1), m);
    call dgemm('N', 'T', m, n  , n-2, 1.0d0, aux1(:,2:n-1), m, V2, n, 1.0d0, uu, m);

    ! Plot the solution
    call surf(m,n,x,y,uu);
end
