%% 2D Finite Difference: Poisson Equation for Fluid Flow in Rectangular Duct
% Problem: Solve ρ = μ(∂²vx/∂y² + ∂²vx/∂z²) for fluid in rectangular duct
% For no pressure gradient: ρ = 0, so ∂²vx/∂y² + ∂²vx/∂z² = 0
% Boundary conditions: vx = 0 (fixed wall) or vx = v0 (moving wall)

clear; clc; close all;

%% Problem Setup
% Domain: rectangular duct [0,Ly] × [0,Lz]
% Equation: ∂²vx/∂y² + ∂²vx/∂z² = 0 (Laplace equation)

%% Grid Parameters
Ny = 41;  % Number of grid points in y-direction
Nz = 41;  % Number of grid points in z-direction
Ly = 1.0; % Duct height
Lz = 1.0; % Duct width
v0 = 1.0; % Velocity of moving walls

% Grid spacing
dy = Ly / (Ny - 1);
dz = Lz / (Nz - 1);
y = linspace(0, Ly, Ny);
z = linspace(0, Lz, Nz);

% Total number of unknowns
N = Ny * Nz;

%% Function to convert 2D indices to 1D
% For node (i,j): index = (j-1)*Ny + i
idx = @(i, j) (j-1)*Ny + i;

%% Matrix Setup for the Laplace Equation
% Discretizing ∂²vx/∂y² + ∂²vx/∂z² = 0
% ∂²vx/∂y² ≈ (vx(i+1,j) - 2*vx(i,j) + vx(i-1,j))/dy²
% ∂²vx/∂z² ≈ (vx(i,j+1) - 2*vx(i,j) + vx(i,j-1))/dz²
% Combined: (vx(i+1,j) - 2*vx(i,j) + vx(i-1,j))/dy² + (vx(i,j+1) - 2*vx(i,j) + vx(i,j-1))/dz² = 0

% Function to set up matrix for given boundary conditions
function [A, b] = setupSystem(Ny, Nz, dy, dz, v_bottom, v_top, v_left, v_right)
    N = Ny * Nz;
    A = sparse(N, N);
    b = zeros(N, 1);
    
    % Convert 2D indices to 1D index
    idx = @(i, j) (j-1)*Ny + i;
    
    for j = 1:Nz
        for i = 1:Ny
            n = idx(i, j);
            
            % Check if this is a boundary node
            if i == 1  % Bottom boundary (y = 0)
                A(n, n) = 1;
                b(n) = v_bottom;
            elseif i == Ny  % Top boundary (y = Ly)
                A(n, n) = 1;
                b(n) = v_top;
            elseif j == 1  % Left boundary (z = 0)
                A(n, n) = 1;
                b(n) = v_left;
            elseif j == Nz  % Right boundary (z = Lz)
                A(n, n) = 1;
                b(n) = v_right;
            else
                % Interior node - apply discretized Laplace equation
                % (vx(i+1,j) - 2*vx(i,j) + vx(i-1,j))/dy² + (vx(i,j+1) - 2*vx(i,j) + vx(i,j-1))/dz² = 0
                
                % Coefficient for vx(i-1,j)
                A(n, idx(i-1, j)) = 1/dy^2;
                
                % Coefficient for vx(i+1,j)
                A(n, idx(i+1, j)) = 1/dy^2;
                
                % Coefficient for vx(i,j-1)
                A(n, idx(i, j-1)) = 1/dz^2;
                
                % Coefficient for vx(i,j+1)
                A(n, idx(i, j+1)) = 1/dz^2;
                
                % Coefficient for vx(i,j)
                A(n, n) = -2 * (1/dy^2 + 1/dz^2);
                
                % Right-hand side is zero for Laplace equation
                b(n) = 0;
            end
        end
    end
end

%% Case 1: One Moving Wall (Top wall moving)
fprintf('Solving Case 1: One Moving Wall (Top)\n');
[A1, b1] = setupSystem(Ny, Nz, dy, dz, 0, v0, 0, 0);  % Bottom=0, Top=v0, Left=0, Right=0
vx1 = A1 \ b1;
vx1_2D = reshape(vx1, Ny, Nz);

%% Case 2: Two Opposite Moving Walls (Top and Bottom)
fprintf('Solving Case 2: Two Opposite Moving Walls (Top and Bottom)\n');
[A2, b2] = setupSystem(Ny, Nz, dy, dz, v0, v0, 0, 0);  % Bottom=v0, Top=v0, Left=0, Right=0
vx2 = A2 \ b2;
vx2_2D = reshape(vx2, Ny, Nz);

%% Case 3: Three Moving Walls (Top, Bottom, and Left)
fprintf('Solving Case 3: Three Moving Walls (Top, Bottom, and Left)\n');
[A3, b3] = setupSystem(Ny, Nz, dy, dz, v0, v0, v0, 0);  % Bottom=v0, Top=v0, Left=v0, Right=0
vx3 = A3 \ b3;
vx3_2D = reshape(vx3, Ny, Nz);

%% Plotting
[Z, Y] = meshgrid(z, y);

% Case 1: One Moving Wall
figure('Position', [100, 100, 1200, 500]);

subplot(1, 2, 1);
contourf(Z, Y, vx1_2D, 20, 'LineColor', 'none');
colorbar;
xlabel('z-coordinate');
ylabel('y-coordinate');
title('Case 1: One Moving Wall (Top) - Contour Plot');
axis equal tight;

subplot(1, 2, 2);
surf(Z, Y, vx1_2D);
colorbar;
xlabel('z-coordinate');
ylabel('y-coordinate');
zlabel('Velocity (vx)');
title('Case 1: One Moving Wall (Top) - Surface Plot');
axis tight;

% Case 2: Two Opposite Moving Walls
figure('Position', [150, 150, 1200, 500]);

subplot(1, 2, 1);
contourf(Z, Y, vx2_2D, 20, 'LineColor', 'none');
colorbar;
xlabel('z-coordinate');
ylabel('y-coordinate');
title('Case 2: Two Opposite Moving Walls - Contour Plot');
axis equal tight;

subplot(1, 2, 2);
surf(Z, Y, vx2_2D);
colorbar;
xlabel('z-coordinate');
ylabel('y-coordinate');
zlabel('Velocity (vx)');
title('Case 2: Two Opposite Moving Walls - Surface Plot');
axis tight;

% Case 3: Three Moving Walls
figure('Position', [200, 200, 1200, 500]);

subplot(1, 2, 1);
contourf(Z, Y, vx3_2D, 20, 'LineColor', 'none');
colorbar;
xlabel('z-coordinate');
ylabel('y-coordinate');
title('Case 3: Three Moving Walls - Contour Plot');
axis equal tight;

subplot(1, 2, 2);
surf(Z, Y, vx3_2D);
colorbar;
xlabel('z-coordinate');
ylabel('y-coordinate');
zlabel('Velocity (vx)');
title('Case 3: Three Moving Walls - Surface Plot');
axis tight;

%% Comparison of all cases
figure('Position', [250, 250, 1200, 400]);

% Contour plots for comparison
subplot(1, 3, 1);
contourf(Z, Y, vx1_2D, 20, 'LineColor', 'none');
colorbar;
xlabel('z-coordinate');
ylabel('y-coordinate');
title('Case 1: One Moving Wall');
axis equal tight;

subplot(1, 3, 2);
contourf(Z, Y, vx2_2D, 20, 'LineColor', 'none');
colorbar;
xlabel('z-coordinate');
ylabel('y-coordinate');
title('Case 2: Two Opposite Moving Walls');
axis equal tight;

subplot(1, 3, 3);
contourf(Z, Y, vx3_2D, 20, 'LineColor', 'none');
colorbar;
xlabel('z-coordinate');
ylabel('y-coordinate');
title('Case 3: Three Moving Walls');
axis equal tight;

sgtitle('Comparison of Velocity Profiles for Different Boundary Conditions');

%% Analysis of results
fprintf('\n=== 2D POISSON EQUATION ANALYSIS ===\n');
fprintf('Grid: %d × %d points\n', Ny, Nz);
fprintf('Domain: [0,%.1f] × [0,%.1f]\n', Ly, Lz);

fprintf('\nMaximum velocities:\n');
fprintf('Case 1 (One moving wall):           vx_max = %.4f\n', max(max(vx1_2D)));
fprintf('Case 2 (Two opposite moving walls): vx_max = %.4f\n', max(max(vx2_2D)));
fprintf('Case 3 (Three moving walls):        vx_max = %.4f\n', max(max(vx3_2D)));

fprintf('\nPhysical Interpretation:\n');
fprintf('- Case 1: Velocity gradient from fixed walls to the moving wall\n');
fprintf('- Case 2: More uniform velocity field due to symmetry\n');
fprintf('- Case 3: Complex flow pattern with three moving boundaries\n');

% Add the verification of the Laplace equation
fprintf('\nVerification (Laplace equation residual at center):\n');
i_center = round(Ny/2); 
j_center = round(Nz/2);

% Function to calculate the Laplacian ∇²vx at a point
laplacian = @(vx_2D, i, j) ((vx_2D(i+1,j) - 2*vx_2D(i,j) + vx_2D(i-1,j))/dy^2) + ...
                           ((vx_2D(i,j+1) - 2*vx_2D(i,j) + vx_2D(i,j-1))/dz^2);

% Check Laplacian at center point for each case
res1 = laplacian(vx1_2D, i_center, j_center);
res2 = laplacian(vx2_2D, i_center, j_center);
res3 = laplacian(vx3_2D, i_center, j_center);

fprintf('Case 1 Laplacian at center: %.2e (should be ≈ 0)\n', res1);
fprintf('Case 2 Laplacian at center: %.2e (should be ≈ 0)\n', res2);
fprintf('Case 3 Laplacian at center: %.2e (should be ≈ 0)\n', res3);