function x = sherman_morrison_solve(L,U,P,u,v,b)

x1 = U \ (L \ (P * b));
x2 = U \ (L \ (P * u));
alpha = 1 - v' * x2;
beta = v' * x1;

x = x1 + (beta / alpha) * x2;

end
