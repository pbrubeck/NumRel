module visual
contains

subroutine plot(n,x,y)
    implicit none
    integer :: n
    integer :: i
    double precision :: x(n), y(n)

    open(UNIT=64,FILE='data/data1d.dat');
    do i=1,n
        write(64,*) x(i), y(i);
    enddo
    close(64);

    call system('gnuplot -p ./data/plot.plt')
end


subroutine surf(m,n,x,y,z)
    implicit none
    integer :: m, n
    integer :: i, j
    double precision :: x(m), y(n), z(m,n)

    open(UNIT=64,FILE='data/data2d.dat');
    do i=1,m
        do j=1,n
            write(64,*) x(i), y(j), z(i,j), z(i,j);
        enddo
        write(64,*);
    enddo
    close(64);

    call system('gnuplot -p ./data/surf.plt')
end


subroutine mesh(m,n,x,y,z)
    implicit none
    integer :: m, n
    integer :: i, j
    double precision :: x(m), y(n), z(m,n)

    open(UNIT=64,FILE='data/data2d.dat');
    do i=1,m
        do j=1,n
            write(64,*) x(i), y(j), z(i,j), z(i,j);
        enddo
        write(64,*);
    enddo
    close(64);

    call system('gnuplot -p ./data/mesh.plt')
end

end module visual
