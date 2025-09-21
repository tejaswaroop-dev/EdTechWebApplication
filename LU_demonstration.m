%% LU Decomposition Demonstration - Circuit Analysis Example
% This script demonstrates the use of myLUDecomposer and mySolver functions
% in solving a practical engineering problem: DC circuit analysis using
% nodal analysis method.
%
% The example shows how LU decomposition can efficiently solve systems
% of linear equations that arise in electrical circuit analysis.

clear; clc; close all;

fprintf('=== LU Decomposition Demonstration: DC Circuit Analysis ===\n\n');

%% Problem Setup: DC Circuit with Multiple Nodes
% Consider a DC circuit with 4 nodes (excluding ground) and various resistors
% We'll use nodal analysis to find the node voltages

% Circuit parameters
R1 = 100;    % Ohms
R2 = 200;    % Ohms  
R3 = 150;    % Ohms
R4 = 300;    % Ohms
R5 = 250;    % Ohms
R6 = 180;    % Ohms

% Voltage sources
Vs1 = 12;    % Volts
Vs2 = 8;     % Volts

% Current sources  
Is1 = 0.05;  % Amperes

fprintf('Circuit Parameters:\n');
fprintf('  Resistors: R1=%.0fΩ, R2=%.0fΩ, R3=%.0fΩ, R4=%.0fΩ, R5=%.0fΩ, R6=%.0fΩ\n', ...
        R1, R2, R3, R4, R5, R6);
fprintf('  Voltage Sources: Vs1=%.0fV, Vs2=%.0fV\n', Vs1, Vs2);
fprintf('  Current Source: Is1=%.3fA\n\n', Is1);

%% Nodal Analysis Matrix Setup
% For a 4-node circuit, we get a 4x4 conductance matrix
% Each equation represents Kirchhoff's Current Law at a node

% Conductances (1/Resistance)
G1 = 1/R1; G2 = 1/R2; G3 = 1/R3; G4 = 1/R4; G5 = 1/R5; G6 = 1/R6;

% Conductance matrix [G] for nodal analysis
% G(i,j) represents the conductance between nodes i and j
G = zeros(4, 4);

% Node 1 equation: (G1+G2+G3)*V1 - G2*V2 - G3*V3 = Vs1*G1 + Is1
G(1,1) = G1 + G2 + G3;
G(1,2) = -G2;
G(1,3) = -G3;
G(1,4) = 0;

% Node 2 equation: -G2*V1 + (G2+G4+G5)*V2 - G5*V4 = 0
G(2,1) = -G2;
G(2,2) = G2 + G4 + G5;
G(2,3) = 0;
G(2,4) = -G5;

% Node 3 equation: -G3*V1 + (G3+G6)*V3 - G6*V4 = Vs2*G6
G(3,1) = -G3;
G(3,2) = 0;
G(3,3) = G3 + G6;
G(3,4) = -G6;

% Node 4 equation: -G5*V2 - G6*V3 + (G5+G6)*V4 = 0
G(4,1) = 0;
G(4,2) = -G5;
G(4,3) = -G6;
G(4,4) = G5 + G6;

% Current vector [I] (right-hand side)
I = [Vs1*G1 + Is1; 0; Vs2*G6; 0];

fprintf('Conductance Matrix G (in Siemens):\n');
fprintf('   Node1    Node2    Node3    Node4\n');
for i = 1:4
    fprintf('Node%d ', i);
    for j = 1:4
        fprintf('%8.4f ', G(i,j));
    end
    fprintf('\n');
end

fprintf('\nCurrent Vector I (in Amperes):\n');
for i = 1:4
    fprintf('Node%d: %8.4f\n', i, I(i));
end

%% Solution using our LU decomposition functions
fprintf('\n--- Solving using our LU decomposition functions ---\n');

% Step 1: LU Decomposition
tic;
[P, L, U] = myLUDecomposer(G);
decomp_time = toc;

% Step 2: Solve system
tic;
V = mySolver(G, I);
solve_time = toc;

fprintf('LU Decomposition time: %.6f seconds\n', decomp_time);
fprintf('System solution time: %.6f seconds\n', solve_time);

%% Display Results
fprintf('\nNode Voltages (using our solver):\n');
for i = 1:4
    fprintf('  V%d = %.4f V\n', i, V(i));
end

% Verify solution
residual = norm(G*V - I);
fprintf('\nSolution verification:\n');
fprintf('  Residual ||GV - I|| = %.2e\n', residual);

%% Comparison with MATLAB's built-in solver
fprintf('\n--- Comparison with MATLAB''s built-in solver ---\n');
tic;
V_matlab = G \ I;
matlab_time = toc;

fprintf('MATLAB solver time: %.6f seconds\n', matlab_time);
fprintf('Speed ratio (ours/MATLAB): %.2f\n', (decomp_time + solve_time)/matlab_time);

% Compare solutions
diff_solutions = norm(V - V_matlab);
fprintf('Solution difference: %.2e\n', diff_solutions);

fprintf('\nNode Voltages (MATLAB solver):\n');
for i = 1:4
    fprintf('  V%d = %.4f V\n', i, V_matlab(i));
end

%% Power Analysis using the solution
fprintf('\n--- Power Analysis ---\n');

% Calculate currents through each resistor using Ohm's law
I_R1 = (Vs1 - V(1)) / R1;
I_R2 = (V(1) - V(2)) / R2;
I_R3 = (V(1) - V(3)) / R3;
I_R4 = V(2) / R4;
I_R5 = (V(2) - V(4)) / R5;
I_R6 = (V(3) - V(4)) / R6;

% Calculate power dissipation in each resistor
P_R1 = I_R1^2 * R1;
P_R2 = I_R2^2 * R2;
P_R3 = I_R3^2 * R3;
P_R4 = I_R4^2 * R4;
P_R5 = I_R5^2 * R5;
P_R6 = I_R6^2 * R6;

total_power = P_R1 + P_R2 + P_R3 + P_R4 + P_R5 + P_R6;

fprintf('Resistor Currents:\n');
fprintf('  I_R1 = %.4f A, I_R2 = %.4f A, I_R3 = %.4f A\n', I_R1, I_R2, I_R3);
fprintf('  I_R4 = %.4f A, I_R5 = %.4f A, I_R6 = %.4f A\n', I_R4, I_R5, I_R6);

fprintf('\nPower Dissipation:\n');
fprintf('  P_R1 = %.4f W, P_R2 = %.4f W, P_R3 = %.4f W\n', P_R1, P_R2, P_R3);
fprintf('  P_R4 = %.4f W, P_R5 = %.4f W, P_R6 = %.4f W\n', P_R4, P_R5, P_R6);
fprintf('  Total Power = %.4f W\n', total_power);

%% Visualization
figure('Position', [100, 100, 1200, 800]);

% Plot 1: Node voltages
subplot(2, 3, 1);
bar(1:4, V, 'FaceColor', [0.2, 0.6, 0.8]);
xlabel('Node Number');
ylabel('Voltage (V)');
title('Node Voltages');
grid on;
for i = 1:4
    text(i, V(i) + 0.1, sprintf('%.3f V', V(i)), ...
         'HorizontalAlignment', 'center', 'FontSize', 10);
end

% Plot 2: Resistor currents
subplot(2, 3, 2);
currents = [I_R1, I_R2, I_R3, I_R4, I_R5, I_R6];
bar(1:6, currents, 'FaceColor', [0.8, 0.4, 0.2]);
xlabel('Resistor Number');
ylabel('Current (A)');
title('Resistor Currents');
grid on;
set(gca, 'XTickLabel', {'R1', 'R2', 'R3', 'R4', 'R5', 'R6'});

% Plot 3: Power dissipation
subplot(2, 3, 3);
powers = [P_R1, P_R2, P_R3, P_R4, P_R5, P_R6];
bar(1:6, powers, 'FaceColor', [0.6, 0.8, 0.3]);
xlabel('Resistor Number');
ylabel('Power (W)');
title('Power Dissipation');
grid on;
set(gca, 'XTickLabel', {'R1', 'R2', 'R3', 'R4', 'R5', 'R6'});

% Plot 4: Conductance matrix visualization
subplot(2, 3, 4);
imagesc(G);
colorbar;
xlabel('Column (Node)');
ylabel('Row (Node)');
title('Conductance Matrix G');
colormap(gca, 'cool');

% Plot 5: LU matrices visualization
subplot(2, 3, 5);
LU_combined = [L, U];
imagesc(LU_combined);
colorbar;
xlabel('Column');
ylabel('Row');
title('L and U Matrices (Side by Side)');
colormap(gca, 'hot');

% Plot 6: Solution comparison
subplot(2, 3, 6);
plot(1:4, V, 'o-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Our Solver');
hold on;
plot(1:4, V_matlab, 's--', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'MATLAB Solver');
xlabel('Node Number');
ylabel('Voltage (V)');
title('Solution Comparison');
legend('Location', 'best');
grid on;

% sgtitle('DC Circuit Analysis using LU Decomposition'); % Not available in Octave
% Add main title manually for Octave compatibility
if exist('sgtitle', 'builtin')
    sgtitle('DC Circuit Analysis using LU Decomposition');
else
    % Octave alternative - use text annotation
    annotation('textbox', [0.4, 0.95, 0.2, 0.05], 'String', ...
               'DC Circuit Analysis using LU Decomposition', ...
               'HorizontalAlignment', 'center', 'FontSize', 14, ...
               'FontWeight', 'bold', 'EdgeColor', 'none');
end

%% Summary
fprintf('\n=== Summary ===\n');
fprintf('The LU decomposition functions successfully solved a 4-node DC circuit:\n');
fprintf('• LU decomposition computed accurately (residual: %.2e)\n', residual);
fprintf('• Solution matches MATLAB''s built-in solver (difference: %.2e)\n', diff_solutions);
fprintf('• Total computation time: %.6f seconds\n', decomp_time + solve_time);
fprintf('• Circuit analysis complete with node voltages and power calculations\n\n');

fprintf('This demonstrates the practical utility of LU decomposition in:\n');
fprintf('• Electrical circuit analysis\n');
fprintf('• Solving large systems of linear equations\n');
fprintf('• Engineering applications requiring numerical stability\n');