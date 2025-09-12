% HW3 Problem 3: 1D Diffusion-Reaction in a Slab
% Solving: d²θ/dη² - Da·θ = 0
% BC: θ(0) = 1, dθ/dη(1) = 0

clear; clc; close all;

%% Parameters
N = 50;                    % Number of grid points
h = 1/N;                   % Grid spacing
eta = linspace(0, 1, N+1); % Grid points
Da_values = [100, 10, 1, 0.1]; % Damköhler numbers
colors = lines(length(Da_values));

%% Main Solution Loop
figure('Position', [100, 100, 1000, 700]);
hold on;

% Storage for error analysis
errors = zeros(length(Da_values), 1);

for idx = 1:length(Da_values)
    Da = Da_values(idx);
    
    fprintf('\nSolving for Da = %.1f\n', Da);
    
    %% Set up the finite difference system
    % We solve for interior points (θ₁ to θₙ₋₁)
    A = zeros(N-1, N-1);
    b = zeros(N-1, 1);
    
    % Interior points: (θᵢ₋₁ - 2θᵢ + θᵢ₊₁)/h² - Da·θᵢ = 0
    % Rearranged: θᵢ₋₁ - (2 + Da·h²)θᵢ + θᵢ₊₁ = 0
    
    for i = 1:N-1
        A(i, i) = -(2 + Da*h^2);  % Diagonal term
        
        if i > 1
            A(i, i-1) = 1;        % Lower diagonal
        end
        
        if i < N-1
            A(i, i+1) = 1;        % Upper diagonal
        end
    end
    
    %% Apply boundary conditions
    
    % Left boundary: θ(0) = 1 (Dirichlet)
    % This affects the first equation: θ₀ - (2 + Da·h²)θ₁ + θ₂ = 0
    % Since θ₀ = 1, we move it to RHS: -(2 + Da·h²)θ₁ + θ₂ = -1
    b(1) = -1;
    
    % Right boundary: dθ/dη(1) = 0 (Neumann)
    % Using second-order backward difference: (3θₙ - 4θₙ₋₁ + θₙ₋₂)/(2h) = 0
    % This gives: 3θₙ - 4θₙ₋₁ + θₙ₋₂ = 0
    A(N-1, N-1) = 3;      % θₙ coefficient
    if N-1 > 1
        A(N-1, N-2) = -4; % θₙ₋₁ coefficient
    end
    if N-1 > 2
        A(N-1, N-3) = 1;  % θₙ₋₂ coefficient
    end
    b(N-1) = 0;           % RHS for Neumann condition
    
    %% Solve the linear system
    theta_interior = A\b;
    
    % Reconstruct full solution including boundary point
    theta_numerical = [1; theta_interior]; % θ₀ = 1
    
    %% Analytical solution
    sqrt_Da = sqrt(Da);
    theta_analytical = cosh(sqrt_Da*(1-eta)) / cosh(sqrt_Da);
    
    %% Calculate error
    error_L2 = sqrt(h * sum((theta_numerical - theta_analytical').^2));
    errors(idx) = error_L2;
    
    fprintf('L2 Error: %.2e\n', error_L2);
    
    %% Plot results
    plot(eta, theta_numerical, '-o', 'Color', colors(idx,:), 'LineWidth', 2, ...
         'MarkerSize', 4, 'MarkerIndices', 1:5:length(eta), ...
         'DisplayName', sprintf('Numerical Da=%.1f', Da));
    plot(eta, theta_analytical, '--', 'Color', colors(idx,:), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('Analytical Da=%.1f', Da));
end

%% Format plot
xlabel('\eta = z/L', 'FontSize', 12);
ylabel('\theta = C_A/C_{A0}', 'FontSize', 12);
title('1D Diffusion-Reaction through a Slab', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 10);
grid on;
xlim([0, 1]);
ylim([0, 1]);
hold off;

%% Display error summary
fprintf('\n=== ERROR ANALYSIS ===\n');
fprintf('Da\t\tL2 Error\n');
fprintf('---\t\t--------\n');
for idx = 1:length(Da_values)
    fprintf('%.1f\t\t%.2e\n', Da_values(idx), errors(idx));
end

%% Grid convergence study
fprintf('\n=== GRID CONVERGENCE STUDY ===\n');
N_values = [25, 50, 100, 200];
Da_test = 1.0;  % Test with Da = 1

figure('Position', [1100, 100, 600, 500]);
errors_conv = zeros(length(N_values), 1);

for i = 1:length(N_values)
    N_test = N_values(i);
    h_test = 1/N_test;
    eta_test = linspace(0, 1, N_test+1);
    
    % Set up system for this grid
    A_test = zeros(N_test-1, N_test-1);
    b_test = zeros(N_test-1, 1);
    
    for j = 1:N_test-1
        A_test(j, j) = -(2 + Da_test*h_test^2);
        if j > 1
            A_test(j, j-1) = 1;
        end
        if j < N_test-1
            A_test(j, j+1) = 1;
        end
    end
    
    b_test(1) = -1;
    A_test(N_test-1, N_test-1) = 3;
    if N_test-1 > 1
        A_test(N_test-1, N_test-2) = -4;
    end
    if N_test-1 > 2
        A_test(N_test-1, N_test-3) = 1;
    end
    b_test(N_test-1) = 0;
    
    % Solve and calculate error
    theta_test = [1; A_test\b_test];
    sqrt_Da_test = sqrt(Da_test);
    theta_exact = cosh(sqrt_Da_test*(1-eta_test)) / cosh(sqrt_Da_test);
    
    errors_conv(i) = sqrt(h_test * sum((theta_test - theta_exact').^2));
    
    fprintf('N = %3d, h = %.4f, Error = %.2e\n', N_test, h_test, errors_conv(i));
end

% Plot convergence
loglog(1./N_values, errors_conv, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
loglog(1./N_values, 0.1*(1./N_values).^2, 'r--', 'LineWidth', 1.5, 'DisplayName', 'h²');
xlabel('Grid spacing h', 'FontSize', 12);
ylabel('L2 Error', 'FontSize', 12);
title(sprintf('Grid Convergence Study (Da = %.1f)', Da_test), 'FontSize', 14);
legend('Numerical Error', 'Second-order slope', 'Location', 'best');
grid on;

%% Physical interpretation
fprintf('\n=== PHYSICAL INTERPRETATION ===\n');
for idx = 1:length(Da_values)
    Da = Da_values(idx);
    penetration_depth = 1/sqrt(Da);
    
    if Da > 10
        regime = 'Reaction-limited';
    elseif Da < 0.1
        regime = 'Diffusion-limited';
    else
        regime = 'Mixed regime';
    end
    
    fprintf('Da = %.1f: %s, Penetration depth = %.3f\n', ...
            Da, regime, penetration_depth);
end

fprintf('\nSolution completed successfully!\n');
