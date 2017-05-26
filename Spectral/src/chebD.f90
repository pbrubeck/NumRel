subroutine chebD(n, D,x)
    implicit none
    integer, intent(in) :: n
    double precision :: x(n), D(n,n)
    double precision pi
    parameter (pi=4*atan(1.0d0))
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
