subroutine poisson(n,alph,beta,bc)
    use visual
    implicit none
    integer, intent(in) :: n
    integer :: i, info, rd(2), ipiv(n-2)
    double precision, intent(in) :: alph(2), beta(2), bc(2)
    double precision :: x(n),u(n),D(n,n),D2(n,n),B(2,n),G(2,n-2),A(n-2,n-2)

    rd=[1,n];
    !kd=[(i, i=2,n-1)];

    ! Differential operators D, D2
    call chebD(n, D, x);
    call dgemm('N', 'N', n, n, n, 1.0d0, D, n, D, n, 0.0d0, D2, n);

    ! Source term
    u=-5*(x*x-0.25d0);

    ! Boundary conditions
    ! B=diag(alph)*E(rd,:)+diag(beta)*D(rd,:);
    B(1,:)=beta(1)*D(1,:);
    B(1,1)=alph(1)+B(1,1);
    B(2,:)=beta(2)*D(n,:);
    B(2,n)=alph(2)+B(2,n);
    ! G=-B(:,rd)\B(:,kd);
    G=B(:,2:n-1);
    call dgesv(2, n-2, -B(:,rd), 2, ipiv(1:2), G, 2, info);
    ! A=D2(kd,kd)+D2(kd,rd)*G;
    A=D2(2:n-1,2:n-1);
    call dgemm('N', 'N', n-2, n-2, 2, 1.0d0, D2(2:n-1,rd), n-2, G, 2, 1.0d0, A, n-2);

    ! Solution
    ! u(rd)=B(:,rd)\bc;
    call dcopy(2, bc, 1, u, rd(2)-rd(1));
    call dgesv(2, 1, B(:,rd), 2, ipiv(1:2), u(1:n:n-1), 2, info);
    ! u(kd)=u(kd)-D2(kd,rd)*u(rd);
    u(2:n-1)=u(2:n-1)-u(1)*D2(2:n-1,1)-u(n)*D2(2:n-1,n);
    ! u(kd)=A\u(kd)
    call dgesv(n-2, 1, A, n-2, ipiv, u(2:n-1), n-2, info);
    ! u(rd)=u(rd)+G*u(kd);
    call dgemv('N', 2, n-2, 1.0d0, G, 2, u(2:n-1), 1, 1.0d0, u, rd(2)-rd(1));

    call plot(n,x,u);
end
