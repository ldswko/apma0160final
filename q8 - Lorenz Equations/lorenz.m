sigma = 10;
r = 28;
b = 8/3;

t_span = [0, 20];
y0 = [5; 5; 5];

lorenz_eqns = @(t, y) [
    sigma * (y(2) - y(1));
    r * y(1) - y(2) - y(1) * y(3);
    y(1) * y(2) - b * y(3)];

[t, y] = ode45(lorenz_eqns, t_span, y0);

x_vals = y(:, 1);
y_vals = y(:, 2);
z_vals = y(:, 3);

plot3(x_vals, y_vals, z_vals, 'LineWidth', 1.5);
grid on;
xlabel('x(t)');
ylabel('y(t)');
zlabel('z(t)');
title('Lorenz Attractor');
view(3);
