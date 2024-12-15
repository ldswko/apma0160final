function dawson_ode

t0 = 0; t_end = 5; 
N1 = 200; 
N2 = 10;   
t_200 = linspace(t0, t_end, N1 + 1);
t_10 = linspace(t0, t_end, N2 + 1);
y0 = 0;
    
y_heun = zeros(size(t_200));
y_rk4 = zeros(size(t_200));
y_trap = zeros(size(t_200));
y_heun_2 = zeros(size(t_10));
y_rk4_2 = zeros(size(t_10));
y_trap_2 = zeros(size(t_10));
    
[y_heun, y_rk4, y_trap] = calculate_dawson(t_200, y0);
[y_heun_2, y_rk4_2, y_trap_2] = calculate_dawson(t_10, y0);
y_exact_200 = dawson_integral(t_200);
y_exact_10 = dawson_integral(t_10);
    
heun_err = max(abs(y_heun - y_exact_200));
rk4_err = max(abs(y_rk4 - y_exact_200));
trap_err = max(abs(y_trap - y_exact_200));
heun_err_2 = max(abs(y_heun_2 - y_exact_10));
rk4_err_2 = max(abs(y_rk4_2 - y_exact_10));
trap_err_2 = max(abs(y_trap_2 - y_exact_10));
    
fprintf('Maximum Errors with 200 Steps:\n');
fprintf('Heun: %f, RK4: %f, Trapezoid: %f\n', heun_err, rk4_err, trap_err);
fprintf('Maximum Errors with 10 Steps:\n');
fprintf('Heun: %f, RK4: %f, Trapezoid: %f\n', heun_err_2, rk4_err_2, trap_err_2);

end

function [y_heun, y_rk4, y_trap] = calculate_dawson(t, y0)

N = length(t) - 1;
y_heun = zeros(size(t));
y_rk4 = zeros(size(t));
y_trap = zeros(size(t));
    
y_heun(1) = y0;
y_rk4(1) = y0;
y_trap(1) = y0;
    
dt = t(2) - t(1);
    
for k = 1:N

    tk = t(k);
    yk_heun = y_heun(k);
    yk_rk4 = y_rk4(k);
    yk_trap = y_trap(k);
        
    f1 = 1 - 2 * tk * yk_heun;
    y_pred = yk_heun + dt * f1; 
    f2 = 1 - 2 * (tk + dt) * y_pred;
    y_heun(k + 1) = yk_heun + (dt / 2) * (f1 + f2); % Corrector
        
    k1 = 1 - 2 * tk * yk_rk4;
    k2 = 1 - 2 * (tk + dt / 2) * (yk_rk4 + (dt / 2) * k1);
    k3 = 1 - 2 * (tk + dt / 2) * (yk_rk4 + (dt / 2) * k2);
    k4 = 1 - 2 * (tk + dt) * (yk_rk4 + dt * k3);
    y_rk4(k + 1) = yk_rk4 + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
        
    y_trap(k + 1) = (yk_trap * (1 - dt * tk) + dt) / (1 + dt * (tk + dt));

end

end
