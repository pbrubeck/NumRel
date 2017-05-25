program main
    implicit none
    call poisson(64,(/1.0d0,1.0d0/),(/0.0d0,0.0d0/),(/1.0d0,-1.0d0/))
    print *, 'Program terminated.'
end program main


subroutine poisson(n,alph,beta,bc)
    implicit none
    integer n
    double precision x(n),u(n),D(n,n),D2(n,n),B(2,n),G(2,n-2),A(n-2,n-2)
    double precision alph(2), beta(2), bc(2)
    integer i, info, rd(2), kd(n-2), ipiv(n-2)
    rd=(/1,n/)

    ! Differential operators D, D2
    call chebD(n,D,x)
    call dgemm('N','N',n,n,n,1.0d0,D,n,D,n,0.0d0,D2,n)

    ! Boundary conditions
    B(1,:)=beta(1)*D(1,:)
    B(1,1)=alph(1)+B(1,1)
    B(2,:)=beta(2)*D(n,:)
    B(2,n)=alph(2)+B(2,n)
    ! G=-B(:,rd)\B(:,kd)
    G=B(:,2:n-1)
    call dgesv(2, n-2, -B(:,rd), 2, ipiv(1:2), G, 2, info)
    ! A=D2(kd,kd)+D2(kd,rd)*G
    A=D2(2:n-1,2:n-1)
    call dgemm('N','N',n-2,n-2,2,1.0d0,D2(2:n-1,rd),n-2,G,2,1.0d0,A,n-2)

    u=-x*x
    u(rd)=bc
    ! u(rd)=-B(:,rd)\bc
    call dgesv(2, 1, -B(:,rd), 2, ipiv(1:2), u(rd), 2, info)
    ! u(kd)=u(kd)-D2(kd,rd)*u(rd)
    u(2:n-1)=u(2:n-1)-u(1)*D2(2:n-1,1)-u(n)*D2(2:n-1,n)
    ! u(kd)=A\u(kd)
    call dgesv(n-2, 1, A, n-2, ipiv, u(2:n-1), n-2, info)
    ! u(rd)=u(rd)+G*u(kd)
    call dgemv('N', 2, n-2, 1.0d0, G, 2, u(2:n-1), 1, 1.0d0, u(rd), 1)

    call disp(u,n,1)
    print *, info
end


subroutine chebD(n, D,x)
    implicit none
    integer, intent(in) :: n
    double precision :: x(n), D(n,n)
    double precision pi
    parameter (pi=4*atan(1.0))
    integer :: i, j, k, ci, cj
    x=(/(cos(pi/(n-1)*i), i=0,n-1)/)
    ci=-1;
    do i=1,n
        ci=sign(merge(2,1,mod(i-1,n-1)==0),-ci);
        cj=-ci;
        do k=i-1,i+n-2
            j=1+mod(k,n)
            cj=sign(merge(2,1,mod(j-1,n-1)==0),-cj);
            D(i,j)=ci/(cj*(x(i)-x(j)+(i==j)));
            D(i,i)=D(i,i)-D(i,j);
        enddo
    enddo
end


subroutine disp(A,m,n)
    implicit none
    integer, intent(in) :: m,n
    double precision, intent(in) :: A(m,n)
    integer :: i,j
    do i=1,size(A,1)
        write(*,"(100f15.5)") ( A(i,j), j=1,size(A,2) )
    enddo
end
