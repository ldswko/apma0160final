function x = forwardsub_diag(l_vals, b)

n = length(b);
x = zeros(n, 1);

x(1) = b(1); 

for i = 2:n
    x(i) = b(i) - l_vals(i-1) * x(i-1);
end

end
