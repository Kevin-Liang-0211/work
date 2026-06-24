clc; clear; close all;

%% Problem 2 (Optional): Helmholtz equation
% phi_xx + phi_yy = lambda^2 phi, (x,y) in (0,1)x(0,1)
% phi(0,y) = sinh(mu*(1-y))/sinh(mu)
% phi(1,y) = -sinh(mu*(1-y))/sinh(mu)
% phi(x,0) = cos(pi*x), phi(x,1) = 0
% mu = sqrt(lambda^2 + pi^2)
% Exact solution: phi(x,y) = cos(pi*x)*sinh(mu*(1-y))/sinh(mu)

h = 0.025;
lambda_list = [0.5, 1.0, 2.0];

fprintf('=== Problem 2 (Optional): Helmholtz equation ===\n');
fprintf('Grid size: h = %.4f\n\n', h);
fprintf('%-10s %-14s\n', 'lambda', 'Max error');
fprintf(repmat('-', 1, 28));
fprintf('\n');

for k = 1:length(lambda_list)
    lambda = lambda_list(k);
    [Phi, x, y, max_err] = helmholtz_five_point_solver(lambda, h);
    fprintf('%-10.2f %-14.6e\n', lambda, max_err);

    figure(k)
    surf(x, y, Phi');
    title(['Helmholtz solution, \lambda=', num2str(lambda), ', h=', num2str(h)]);
    xlabel('x'); ylabel('y'); zlabel('\phi');
    set(gca, 'FontSize', 12);
    view(-55, 25);
    print('-dpng', '-r300', ['CH4_HW2_lambda_', strrep(num2str(lambda), '.', 'p'), '.png']);
end
