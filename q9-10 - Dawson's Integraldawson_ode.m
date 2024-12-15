function dawson_ode

t_full = linspace(0, 5, 201);
t_reduced = linspace(0, 5, 11);
y0 = 0;
f = @(t, y) 1 - 2 * t .* y;

y_heun = heun_method(f, t_full, y0);
y_rk4 = rk4_method(f, t_full, y0);
y_trap = implicit_trapezoid(f, t_full, y0);

y_heun_2 = heun_method(f, t_reduced, y0);
y_rk4_2 = rk4_method(f, t_reduced, y0);
y_trap_2 = implicit_trapezoid(f, t_reduced, y0);

exact_full = dawson_integral(t_full);
exact_reduced = dawson_integral(t_reduced);

heun_err = max(abs(y_heun - exact_full));
rk4_err = max(abs(y_rk4 - exact_full));
trap_err = max(abs(y_trap - exact_full));

heun_err_2 = max(abs(y_heun_2 - exact_reduced));
rk4_err_2 = max(abs(y_rk4_2 - exact_reduced));
trap_err_2 = max(abs(y_trap_2 - exact_reduced));

fprintf('Full time steps over 200 intervals:\n');
fprintf('Heun error: %.4f\n', heun_err);
fprintf('RK4 error: %.4f\n', rk4_err);
fprintf('Trapezoid error: %.4f\n\n', trap_err);
    
fprintf('Reduced time steps over 10 intervals:\n');
fprintf('Heun error: %.4f\n', heun_err_2);
fprintf('RK4 error: %.4f\n', rk4_err_2);
fprintf('Trapezoid error: %.4f\n', trap_err_2);

end

function y = heun_method(f, t, y0)

dt = t(2) - t(1);
y = zeros(size(t));
y(1) = y0;

for k = 1:(length(t)-1)
    y_predict = y(k) + dt * f(t(k), y(k));
    y(k+1) = y(k) + (dt / 2) * (f(t(k), y(k)) + f(t(k+1), y_predict));
end

y = y(:)';

end

function y = rk4_method(f, t, y0)

dt = t(2) - t(1);
y = zeros(size(t));
y(1) = y0;

for k = 1:(length(t)-1)
    k1 = f(t(k), y(k));
    k2 = f(t(k) + dt/2, y(k) + dt/2 * k1);
    k3 = f(t(k) + dt/2, y(k) + dt/2 * k2);
    k4 = f(t(k) + dt, y(k) + dt * k3);
    y(k+1) = y(k) + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4);
end

y = y(:)';

end

function y = implicit_trapezoid(f, t, y0)

dt = t(2) - t(1);
y = zeros(size(t));
y(1) = y0;

for k = 1:(length(t)-1)
    y_next = y(k);
    for iter = 1:10
        y_next = y(k) + (dt / 2) * (f(t(k), y(k)) + f(t(k+1), y_next));
    end
    y(k+1) = y_next;
end

y = y(:)';
 
end
