subroutine chebD(n, D,x)
    implicit none
    integer, intent(in) :: n
    integer :: i, j, c(n)
    double precision :: x(n), D(n,n)
    double precision, parameter :: pi=4*atan(1.0d0)
    x=[(cos(pi/(n-1)*i), i=0,n-1)];
    c=[(sign(merge(2,1,mod(i,n-1)==0),1-2*mod(i,2)), i=0,n-1)];
    do i=1,n
        D(i,:)=[(c(i)/(c(j)*(x(i)-x(j)+(i==j))), j=1,n)];
        D(i,i)=D(i,i)-sum(D(i,:));
    enddo
end
