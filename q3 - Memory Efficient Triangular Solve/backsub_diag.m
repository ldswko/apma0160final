function x = backsub_diag(u_diag, u_super, b)

n = length(b); 
x = zeros(n, 1);

x(n) = b(n) / u_diag(n); 

for i = n-1:-1:1
    x(i) = (b(i) - u_super(i) * x(i+1)) / u_diag(i);
end

end
