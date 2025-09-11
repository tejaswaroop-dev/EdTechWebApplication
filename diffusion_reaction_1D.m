%% 1D Finite Differences: Species Diffusion with First-Order Reaction
% Problem: Species A diffuses through a slab with consumption reaction
% Governing equation: 0 = D * d²CA/dz² - k*CA
% Boundary conditions: CA(0) = CA0, dCA/dz|z=L = 0

clear; clc; close all;

%% Problem Setup and Nondimensionalization
% Original equation: 0 = D * d²CA/dz² - k*CA
% Let: ζ = z/L (dimensionless position)
%      C = CA/CA0 (dimensionless concentration)
% 
% Substituting: d²CA/dz² = (CA0/L²) * d²C/dζ²
% 
% Nondimensional equation: 0 = D*(CA0/L²)*d²C/dζ² - k*CA0*C
% Dividing by k*CA0: 0 = (D/(k*L²))*d²C/dζ² - C
% 
% Define Damkohler number: Da = k*L²/D
% Final nondimensional equation: d²C/dζ² = Da*C
% 
% Boundary conditions:
% C(0) = 1 (at ζ = 0)
% dC/dζ|ζ=1 = 0 (at ζ = 1)

%% Numerical Parameters
N = 50;  % Number of grid points
Da_values = [100, 10, 1, 0.1];  % Damkohler numbers
colors = ['r', 'b', 'g', 'm'];  % Colors for each Da value

% Grid setup
zeta = linspace(0, 1, N);  % Dimensionless position
h = zeta(2) - zeta(1);     % Grid spacing

%% Analytical Solution
% For d²C/dζ² = Da*C with C(0) = 1 and dC/dζ|ζ=1 = 0
% General solution: C(ζ) = A*exp(sqrt(Da)*ζ) + B*exp(-sqrt(Da)*ζ)
% 
% Applying BC1: C(0) = 1 → A + B = 1
% Applying BC2: dC/dζ|ζ=1 = 0 → sqrt(Da)*[A*exp(sqrt(Da)) - B*exp(-sqrt(Da))] = 0
% This gives: A*exp(sqrt(Da)) = B*exp(-sqrt(Da))
% 
% Solving: B = A*exp(2*sqrt(Da))
% From A + B = 1: A = 1/(1 + exp(2*sqrt(Da)))
% Therefore: B = exp(2*sqrt(Da))/(1 + exp(2*sqrt(Da)))

analytical_solution = @(zeta, Da) ...
    (exp(sqrt(Da)*zeta) + exp(sqrt(Da)*(2-zeta))) / (1 + exp(2*sqrt(Da)));

%% Setup and Solve Numerical System for Each Da
figure('Position', [100, 100, 1000, 600]);

for idx = 1:length(Da_values)
    Da = Da_values(idx);
    
    %% Matrix Setup for Finite Differences
    % Interior points: d²C/dζ² ≈ (C_{i+1} - 2*C_i + C_{i-1})/h²
    % Equation becomes: (C_{i+1} - 2*C_i + C_{i-1})/h² = Da*C_i
    % Rearranging: C_{i+1} - (2 + Da*h²)*C_i + C_{i-1} = 0
    
    A = zeros(N, N);
    b = zeros(N, 1);
    
    % First node (boundary condition): C(0) = 1
    A(1, 1) = 1;
    b(1) = 1;
    
    % Interior nodes (i = 2 to N-1)
    for i = 2:N-1
        A(i, i-1) = 1;                    % C_{i-1}
        A(i, i) = -(2 + Da*h^2);         % C_i
        A(i, i+1) = 1;                   % C_{i+1}
        b(i) = 0;
    end
    
    % Last node (Neumann BC): dC/dζ = 0 at ζ = 1
    % Using backward difference: (C_N - C_{N-1})/h = 0 → C_N = C_{N-1}
    A(N, N-1) = -1;
    A(N, N) = 1;
    b(N) = 0;
    
    %% Solve Linear System
    C_numerical = A \ b;
    
    %% Analytical Solution
    C_analytical = analytical_solution(zeta, Da);
    
    %% Plotting
    subplot(2, 2, idx);
    hold on;
    
    % Numerical solution (solid line with square markers)
    plot(zeta, C_numerical, 's-', 'Color', colors(idx), 'LineWidth', 2, ...
         'MarkerSize', 6, 'MarkerFaceColor', colors(idx), 'DisplayName', 'Numerical');
    
    % Analytical solution (dashed line with circle markers)
    plot(zeta, C_analytical, 'o--', 'Color', colors(idx), 'LineWidth', 2, ...
         'MarkerSize', 6, 'MarkerFaceColor', 'none', 'DisplayName', 'Analytical');
    
    % Formatting
    xlabel('Dimensionless Position (ζ = z/L)');
    ylabel('Dimensionless Concentration (C = C_A/C_{A0})');
    title(sprintf('Da = %.1f', Da));
    legend('Location', 'best');
    grid on;
    
    % Calculate and display maximum error
    max_error = max(abs(C_numerical - C_analytical'));
    text(0.05, 0.95, sprintf('Max Error: %.2e', max_error), ...
         'Units', 'normalized', 'BackgroundColor', 'white');
end

sgtitle('1D Diffusion-Reaction: Numerical vs Analytical Solutions');

%% Summary Plot - All Da values together
figure('Position', [150, 150, 800, 600]);
hold on;

for idx = 1:length(Da_values)
    Da = Da_values(idx);
    
    % Solve numerical system
    A = zeros(N, N);
    b = zeros(N, 1);
    
    A(1, 1) = 1;
    b(1) = 1;
    
    for i = 2:N-1
        A(i, i-1) = 1;
        A(i, i) = -(2 + Da*h^2);
        A(i, i+1) = 1;
        b(i) = 0;
    end
    
    A(N, N-1) = -1;
    A(N, N) = 1;
    b(N) = 0;
    
    C_numerical = A \ b;
    C_analytical = analytical_solution(zeta, Da);
    
    % Plot with same color for matching Da
    plot(zeta, C_numerical, 's-', 'Color', colors(idx), 'LineWidth', 2, ...
         'MarkerSize', 4, 'DisplayName', sprintf('Numerical Da=%.1f', Da));
    plot(zeta, C_analytical, 'o--', 'Color', colors(idx), 'LineWidth', 2, ...
         'MarkerSize', 4, 'MarkerFaceColor', 'none', ...
         'DisplayName', sprintf('Analytical Da=%.1f', Da));
end

xlabel('Dimensionless Position (ζ = z/L)');
ylabel('Dimensionless Concentration (C = C_A/C_{A0})');
title('1D Diffusion-Reaction: Effect of Damkohler Number');
legend('Location', 'best');
grid on;

%% Display Results Summary
fprintf('\n=== 1D DIFFUSION-REACTION ANALYSIS ===\n');
fprintf('Grid points: %d\n', N);
fprintf('Grid spacing: %.4f\n', h);
fprintf('\nDamkohler Number Effects:\n');
for idx = 1:length(Da_values)
    Da = Da_values(idx);
    C_end_analytical = analytical_solution(1, Da);
    fprintf('Da = %6.1f: C(ζ=1) = %.4f\n', Da, C_end_analytical);
end

fprintf('\nPhysical Interpretation:\n');
fprintf('- High Da: Fast reaction → steep concentration gradients\n');
fprintf('- Low Da:  Slow reaction → nearly uniform concentration\n');
fprintf('- Da = kL²/D represents reaction time scale vs diffusion time scale\n');
