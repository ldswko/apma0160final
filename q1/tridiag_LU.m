function [L, U] = tridiag_LU(A)

n = size(A, 1);
L = eye(n);
U = zeros(n); 

U(1, 1) = A(1, 1);
U(1, 2) = A(1, 2);

for i = 2:n
    L(i, i-1) = A(i, i-1) / U(i-1, i-1);
    U(i, i) = A(i, i) - L(i, i-1) * U(i-1, i);
    if i < n
        U(i, i+1) = A(i, i+1);
    end
end

end
