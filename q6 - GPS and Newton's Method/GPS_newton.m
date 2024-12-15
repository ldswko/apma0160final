S = [15600, 7540, 20140, 0.07074;
     18760, 2750, 18610, 0.0722;
     17610, 14630, 13480, 0.0769;
     19170, 610, 18390, 0.07242];
c = 299792.458; 

f = @(X) f_S(X, S);
Jf = @(X) Jf_S(X, S);
X_initial = [0; 0; -6370; 0];
X_new = X_initial;
X_hist = X_initial;

tol = 1e-7; 
max_iter = 10;

for k = 1:max_iter
    
    F = f(X_new);
    J = Jf(X_new);
    
    delta_X = -J \ F;
    X_new = X_new + delta_X;
    X_hist = [X_hist, X_new];
    
    if norm(delta_X, 2) < tol
        break;
    end
    
end

plot_GPS(S, X_hist(:, end));

disp('Estimated location and time delay:');
disp(X_hist(:, end));
