load price_data.mat

t_plot = linspace(times(1), times(end), 1000);

n = length(times);
A_1 = ones(n,2);
A_2 = ones(n,3);
A_5 = ones(n,6);

for i = 1:5
    A_5(:, i) = times.^(6-i);
end

for j = 1:2
    A_2(:, j) = times.^(3-j);
end

for k = 1:1
    A_1(:, k) = times.^(2-k);
end

coeffs_1 = A_1\prices;
coeffs_2 = A_2\prices;
coeffs_5 = A_5\prices;

best_time = -0.5 * (coeffs_2(2) / coeffs_2(1));
best_price = polyval(coeffs_2, best_time);

plotA1 = polyval(coeffs_1, t_plot);
plotA2 = polyval(coeffs_2, t_plot);
plotA5 = polyval(coeffs_5, t_plot);

plot(times, prices, 'k.', 'MarkerSize', 20) 
hold on;
plot(best_time, best_price, 'm.', 'MarkerSize', 20); 
hold on;
plot(t_plot, plotA1, 'b', 'LineWidth', 1); 
hold on;
plot(t_plot, plotA2, 'r', 'LineWidth', 1); 
hold on;
plot(t_plot, plotA5, 'g', 'LineWidth', 1); 

legend('Observed Prices', 'Optimal Time/Price', 'Linear Fit (Degree 1)', ...
    'Quadratic Fit (Degree 2)', 'Polynomial Fit (Degree 5)', ...
    'FontSize', 12, 'Location', 'Best');

xlabel('Times (Days)', 'FontSize', 15);
ylabel('Prices (Dollars)', 'FontSize', 15);

hold off;
