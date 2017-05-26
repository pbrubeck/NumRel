subroutine disp(A,m,n)
    implicit none
    integer, intent(in) :: m,n
    double precision, intent(in) :: A(m,n)
    integer :: i,j
    do i=1,size(A,1)
        write(*,"(100g15.5)") ( A(i,j), j=1,size(A,2) )
    enddo
end
