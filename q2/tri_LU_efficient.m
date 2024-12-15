function [l_sub, u_diag, u_super] = tri_LU_efficient(a_diag, a_sub, a_super)

n = length(a_diag);
l_sub = zeros(1, n-1);
u_diag = zeros(1, n);
u_super = zeros(1, n-1);

u_diag(1) = a_diag(1); 
u_super(1) = a_super(1); 

for i = 1:n-1
    l_sub(i) = a_sub(i) / u_diag(i);
    u_diag(i+1) = a_diag(i+1) - l_sub(i) * a_super(i);
    if i< n-1
        u_super(i+1) = a_super(i+1);
    end
end

end
