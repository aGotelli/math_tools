close all
clear
clc


analytic_function = load("analytic_functions.csv");
interpolated_function = load("interpolated_function.csv");
observed_function = load("observed_function.csv");

Chebyshev_points = load("Chebyshev_points.csv");
interpolation_points = load("interpolation_points.csv");


figure("Name", "f=x")
plot(interpolation_points, analytic_function(1, :), '-r', 'DisplayName', 'analytic');
hold on
plot(interpolation_points, interpolated_function(1, :), '--b', 'DisplayName', 'chebyshev');
plot(Chebyshev_points, observed_function(1, :), 'kx', 'DisplayName', 'chebyshev observations');
legend("Location","northwest")


figure("Name", "f=x²")
plot(interpolation_points, analytic_function(2, :), '-r', 'DisplayName', 'analytic');
hold on
plot(interpolation_points, interpolated_function(2, :), '--b', 'DisplayName', 'chebyshev');
plot(Chebyshev_points, observed_function(2, :), 'kx', 'DisplayName', 'chebyshev observations');
legend("Location","northwest")


figure("Name", "f=x² + sin(x)")
plot(interpolation_points, analytic_function(3, :), '-r', 'DisplayName', 'analytic');
hold on
plot(interpolation_points, interpolated_function(3, :), '--b', 'DisplayName', 'chebyshev');
plot(Chebyshev_points, observed_function(3, :), 'kx', 'DisplayName', 'chebyshev observations');
legend("Location","northwest")




figure("Name", "f=e^x")
plot(interpolation_points, analytic_function(4, :), '-r', 'DisplayName', 'analytic');
hold on
plot(interpolation_points, interpolated_function(4, :), '--b', 'DisplayName', 'chebyshev');
plot(Chebyshev_points, observed_function(4, :), 'kx', 'DisplayName', 'chebyshev observations');
legend("Location","northwest")


figure("Name", "f=sqrt(x)")
plot(interpolation_points, analytic_function(5, :), '-r', 'DisplayName', 'analytic');
hold on
plot(interpolation_points, interpolated_function(5, :), '--b', 'DisplayName', 'chebyshev');
plot(Chebyshev_points, observed_function(5, :), 'kx', 'DisplayName', 'chebyshev observations');
legend("Location","northwest")


figure("Name", "f=atan(x)")
plot(interpolation_points, analytic_function(6, :), '-r', 'DisplayName', 'analytic');
hold on
plot(interpolation_points, interpolated_function(6, :), '--b', 'DisplayName', 'chebyshev');
plot(Chebyshev_points, observed_function(6, :), 'kx', 'DisplayName', 'chebyshev observations');
legend("Location","northwest")



figure("Name", "f= gaussian(x)")
plot(interpolation_points, analytic_function(7, :), '-r', 'DisplayName', 'analytic');
hold on
plot(interpolation_points, interpolated_function(7, :), '--b', 'DisplayName', 'chebyshev');
plot(Chebyshev_points, observed_function(7, :), 'kx', 'DisplayName', 'chebyshev observations');
legend("Location","northwest")
