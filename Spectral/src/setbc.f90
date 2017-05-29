subroutine setbc(n,alph,beta,D,A, V,L,B,H)
    ! Sets boundary conditions for a second-order differential operator.
    ! alph*u + beta*u' = g
    !
    ! Input args
    ! n        : number of collocation points
    ! alph(2)  : boundary condition coefficients for u
    ! beta(2)  : boundary condition coefficients for u'
    ! D(n,n)   : first-order differentiation matrix
    ! A(n,n)   : second-order differential operator
    !
    ! Output args
    ! V(n,n-2)       : homogenous eigenfunctions
    ! L(n-2)         : corresponding eigenvalues
    ! B(2,n)         : boundary condition operator
    ! H(n,2)         : boundary condition dual basis

    use linalg
    implicit none
    integer, intent(in) :: n
    integer :: i, lwork, info, rd(2), ipiv(2)
    double precision, intent(in) :: alph(2), beta(2), D(n,n), A(n,n)
    double precision :: B(2,n), G(2,n-2), S(n-2,n-2)
    double precision :: V(n,n-2), L(n), H(n,2)
    double precision :: LI(n), WORK(8*n), DUMMY(1,1)

    lwork=8*n;
    rd=[1,n];
    !kd=[(i, i=2,n-1)];

    ! Boundary condition operator
    ! B=diag(alph)*E(rd,:)+diag(beta)*D(rd,:);
    B(1,:)=beta(1)*D(1,:);
    B(1,1)=alph(1)+B(1,1);
    B(2,:)=beta(2)*D(n,:);
    B(2,n)=alph(2)+B(2,n);

    ! Boundary condition basis (dual B*H=I)
    ! H(rd,:)=inv(B(:,rd));
    H(rd,:)=matinv2(B(:,rd));

    ! Give-back matrix
    ! G=-B(:,rd)\B(:,kd);
    G=B(:,2:n-1);
    call dgesv(2, n-2, -B(:,rd), 2, ipiv, G, 2, info);

    ! Schur complement (Sigma, A tilde)
    ! S=A(kd,kd)+A(kd,rd)*G;
    S=A(2:n-1,2:n-1);
    call dgemm('N', 'N', n-2, n-2, 2, 1.0d0, A(2:n-1,1:n:n-1), n-2, G, 2, 1.0d0, S, n-2);

    ! Eigenvalues and eigenvectors
    ! [V(kd,:),L(kd)]=eig(S,'vector');
    call dgeev('N', 'V', n-2, S, n-2, L(2:n-1), LI(2:n-1), DUMMY, 1, V(2:n-1,:), n-2, WORK, lwork, info);
    ! V(rd,:)=G*V(kd,:);
    call dgemm('N', 'N', 2, n-2, n-2, 1.0d0, G, 2, V(2:n-1,:), n-2, 0.0d0, V(1:n:n-1,:), 2);
end
