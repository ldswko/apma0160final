function root = newton_bisection(a, b, f, df, tol, max_iter)

fa = f(a);
fb = f(b);

if fa * fb >= 0
    error('No sign change on interval.');
end

x_old = a;

for k = 1:max_iter

    x_newton = x_old - f(x_old) / df(x_old); 
    if a < x_newton && x_newton < b
        x_new = x_newton;
    else
        x_new = (a + b) / 2;
    end
    if abs(x_new - x_old) < tol
        root = x_new;
        return;
    end

    if f(a) * f(x_new) < 0
        b = x_new;
    else
        a = x_new;
    end
    
    x_old = x_new;

end
    
root = x_new;

end
