function out = Jf_S(X, S)

c = 299792.458; 
x = X(1); 
y = X(2); 
z = X(3); 
d = X(4);

out = zeros(4, 4); 

for i = 1:4
    Ai = S(i, 1);
    Bi = S(i, 2);
    Ci = S(i, 3);
    
    denom = sqrt((x - Ai)^2 + (y - Bi)^2 + (z - Ci)^2);
    
    out(i, 1) = (x - Ai) / denom; 
    out(i, 2) = (y - Bi) / denom; 
    out(i, 3) = (z - Ci) / denom; 
    out(i, 4) = c;  
    
end

end
