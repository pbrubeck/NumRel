subroutine helmholtz(n,alph,beta, V,L)
    implicit none
    integer, intent(in) :: n
    integer :: i, info, rd(2),  ipiv(n-2)
    double precision, intent(in) :: alph(2), beta(2)
    double precision :: x(n),D(n,n),D2(n,n),B(2,n),G(2,n-2),A(n-2,n-2)
    double precision :: V(n,n-2),L(n),LI(n),WORK(8*n),DUMMY(1,1)

    rd=[1,n]
    !kd=(/(i, i=2,n-1)/)

    ! Differential operators D, D2
    call chebD(n, D, x)
    call dgemm('N', 'N', n, n, n, 1.0d0, D, n, D, n, 0.0d0, D2, n)

    ! Boundary conditions
    B(1,:)=beta(1)*D(1,:)
    B(1,1)=alph(1)+B(1,1)
    B(2,:)=beta(2)*D(n,:)
    B(2,n)=alph(2)+B(2,n)

    ! G=-B(:,rd)\B(:,kd);
    G=B(:,2:n-1)
    call dgesv(2, n-2, -B(:,rd), 2, ipiv(1:2), G, 2, info)
    ! A=D2(kd,kd)+D2(kd,rd)*G;
    A=D2(2:n-1,2:n-1)
    call dgemm('N', 'N', n-2, n-2, 2, 1.0d0, D2(2:n-1,1:n:n-1), n-2, G, 2, 1.0d0, A, n-2)
    ! [V(kd,:),L(kd)]=eig(A,'vector');
    call dgeev('N', 'V', n-2, A, n-2, L(2:n-1), LI(2:n-1), DUMMY, 1, V(2:n-1,:), n-2, WORK, 8*n, info)
    ! V(rd,:)=G*V(kd,:);
    call dgemm('N', 'N', 2, n-2, n-2, 1.0d0, G, 2, V(2:n-1,:), n-2, 0.0d0, V(1:n:n-1,:), 2)
end
