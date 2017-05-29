module linalg
contains

pure function matinv2(A) result(B)
    !! Performs a direct calculation of the inverse of a 2×2 matrix.
    double precision, intent(in) :: A(2,2)   !! Matrix
    double precision             :: B(2,2)   !! Inverse matrix
    double precision             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * A(2,2)
    B(2,1) = -detinv * A(2,1)
    B(1,2) = -detinv * A(1,2)
    B(2,2) = +detinv * A(1,1)
end function

pure function sylvester2(A,B,C) result(X)
    !! Solves Sylvester's equation for 2×2 matrices
    double precision, intent(in) :: A(2,2), B(2,2), C(2,2)
    double precision :: detB, trB, cofB(2,2), E(2,2), X(2,2)
    detB=B(1,1)*B(2,2)-B(1,2)*B(2,1);
    trB=B(1,1)+B(2,2);
    cofB(1,1)= B(2,2);
    cofB(2,1)=-B(2,1);
    cofB(1,2)=-B(1,2);
    cofB(2,2)= B(1,1);
    E(:,:)=0.0d0;
    E(1,1)=detB;
    E(2,2)=detB;
    X=matmul(matinv2(E+trB*A+matmul(A, A)), matmul(A, C)+matmul(C, cofB));
end function
end module linalg
