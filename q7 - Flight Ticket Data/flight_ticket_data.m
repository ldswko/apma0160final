load price_data.mat

t_plot = linspace(times(1), times(end), 1000);
A1 = [times, ones(size(times))];
A2 = [times.^2, times, ones(size(times))]; 
A5 = [times.^5, times.^4, times.^3, times.^2, times, ones(size(times))]; 

coeffs_1 = A1 \ prices; 
coeffs_2 = A2 \ prices;
coeffs_5 = A5 \ prices; 

c1 = coeffs_2(1);
c2 = coeffs_2(2);
c3 = coeffs_2(3);

best_time = -c2 / (2 * c1);

best_price = c1 * best_time^2 + c2 * best_time + c3;

p1 = coeffs_1(1) + coeffs_1(2) * t_plot; 
p2 = coeffs_2(1) * t_plot.^2 + coeffs_2(2) * t_plot + coeffs_2(3); 
p5 = coeffs_5(1) * t_plot.^5 + coeffs_5(2) * t_plot.^4 + ...
     coeffs_5(3) * t_plot.^3 + coeffs_5(4) * t_plot.^2 + ...
     coeffs_5(5) * t_plot + coeffs_5(6); 

figure;
plot(times, prices, 'k.', 'MarkerSize', 15); 
hold on;
plot(best_time, best_price, 'ro', 'MarkerSize', 10, 'LineWidth', 2); 
plot(t_plot, p1, 'b--', 'LineWidth', 1.5); 
plot(t_plot, p2, 'g-', 'LineWidth', 2); 
plot(t_plot, p5, 'm-.', 'LineWidth', 1.5); 

legend('Price Data (Observed)', 'Optimal Purchase Point', 'Linear Fit (Degree 1)', ...
       'Quadratic Fit (Degree 2)', 'Polynomial Fit (Degree 5)', 'Location', 'best');
xlabel('Time (days before flight)');
ylabel('Price ($)');
title('Ticket Prices vs. Time with Polynomial Approximations');
grid on; 
hold off;

save('flight_ticket.data.m', 'A1', 'A2', 'A5', 'coeffs_1', 'coeffs_2', 'coeffs_5', 'best_time', 'best_price');
