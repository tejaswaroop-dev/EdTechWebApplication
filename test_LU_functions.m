%% Test Script for myLUDecomposer and mySolver
% This script tests the LU decomposition and linear solver functions
% with various test cases to ensure correctness and numerical stability.

clear; clc; close all;

fprintf('=== Testing myLUDecomposer and mySolver ===\n\n');

%% Test 1: Simple 3x3 matrix
fprintf('Test 1: Simple 3x3 matrix\n');
A1 = [2, 1, 1; 4, 3, 3; 8, 7, 9];
b1 = [4; 10; 24];

try
    % Test LU decomposition
    [P1, L1, U1] = myLUDecomposer(A1);
    
    % Verify P*A = L*U
    error1 = norm(P1*A1 - L1*U1, 'fro');
    fprintf('  LU Decomposition error: %.2e\n', error1);
    
    % Test solver
    x1 = mySolver(A1, b1);
    
    % Verify solution
    residual1 = norm(A1*x1 - b1);
    fprintf('  Solution residual: %.2e\n', residual1);
    fprintf('  Solution: x = [%.4f, %.4f, %.4f]\n', x1(1), x1(2), x1(3));
    
    % Compare with MATLAB's built-in solver
    x1_matlab = A1 \ b1;
    diff1 = norm(x1 - x1_matlab);
    fprintf('  Difference from MATLAB solver: %.2e\n', diff1);
    
catch ME
    fprintf('  ERROR: %s\n', ME.message);
end

fprintf('\n');

%% Test 2: Well-conditioned 4x4 matrix
fprintf('Test 2: Well-conditioned 4x4 matrix\n');
A2 = [4, 3, 2, 1; 3, 4, 3, 2; 2, 3, 4, 3; 1, 2, 3, 4];
b2 = [10; 12; 12; 10];

try
    % Test LU decomposition
    [P2, L2, U2] = myLUDecomposer(A2);
    
    % Verify P*A = L*U
    error2 = norm(P2*A2 - L2*U2, 'fro');
    fprintf('  LU Decomposition error: %.2e\n', error2);
    
    % Test solver
    x2 = mySolver(A2, b2);
    
    % Verify solution
    residual2 = norm(A2*x2 - b2);
    fprintf('  Solution residual: %.2e\n', residual2);
    
    % Compare with MATLAB's built-in solver
    x2_matlab = A2 \ b2;
    diff2 = norm(x2 - x2_matlab);
    fprintf('  Difference from MATLAB solver: %.2e\n', diff2);
    
catch ME
    fprintf('  ERROR: %s\n', ME.message);
end

fprintf('\n');

%% Test 3: Matrix requiring pivoting
fprintf('Test 3: Matrix requiring pivoting\n');
A3 = [0, 1, 2; 1, 2, 3; 2, 3, 1];
b3 = [5; 6; 5];

try
    % Test LU decomposition
    [P3, L3, U3] = myLUDecomposer(A3);
    
    % Verify P*A = L*U
    error3 = norm(P3*A3 - L3*U3, 'fro');
    fprintf('  LU Decomposition error: %.2e\n', error3);
    
    % Test solver
    x3 = mySolver(A3, b3);
    
    % Verify solution
    residual3 = norm(A3*x3 - b3);
    fprintf('  Solution residual: %.2e\n', residual3);
    
    % Compare with MATLAB's built-in solver
    x3_matlab = A3 \ b3;
    diff3 = norm(x3 - x3_matlab);
    fprintf('  Difference from MATLAB solver: %.2e\n', diff3);
    
catch ME
    fprintf('  ERROR: %s\n', ME.message);
end

fprintf('\n');

%% Test 4: Multiple right-hand sides
fprintf('Test 4: Multiple right-hand sides\n');
A4 = [3, 2, 1; 1, 3, 2; 2, 1, 3];
B4 = [6, 14; 6, 11; 6, 13];  % Two RHS vectors

try
    % Test solver with multiple RHS
    X4 = mySolver(A4, B4);
    
    % Verify solutions
    residual4 = norm(A4*X4 - B4, 'fro');
    fprintf('  Solution residual (Frobenius): %.2e\n', residual4);
    
    % Compare with MATLAB's built-in solver
    X4_matlab = A4 \ B4;
    diff4 = norm(X4 - X4_matlab, 'fro');
    fprintf('  Difference from MATLAB solver: %.2e\n', diff4);
    
catch ME
    fprintf('  ERROR: %s\n', ME.message);
end

fprintf('\n');

%% Test 5: Error handling - Singular matrix
fprintf('Test 5: Error handling - Singular matrix\n');
A5 = [1, 2, 3; 2, 4, 6; 1, 1, 1];  % Singular matrix
b5 = [1; 2; 1];

try
    [P5, L5, U5] = myLUDecomposer(A5);
    fprintf('  ERROR: Should have detected singular matrix!\n');
catch ME
    fprintf('  Correctly detected singular matrix: %s\n', ME.message);
end

fprintf('\n');

%% Test 6: Error handling - Non-square matrix
fprintf('Test 6: Error handling - Non-square matrix\n');
A6 = [1, 2, 3; 4, 5, 6];  % 2x3 matrix
b6 = [1; 2];

try
    [P6, L6, U6] = myLUDecomposer(A6);
    fprintf('  ERROR: Should have detected non-square matrix!\n');
catch ME
    fprintf('  Correctly detected non-square matrix: %s\n', ME.message);
end

fprintf('\n');

%% Test 7: Performance comparison with MATLAB's built-in
fprintf('Test 7: Performance comparison (100x100 random matrix)\n');
rng(42);  % For reproducible results
A7 = randn(100, 100) + 10*eye(100);  % Well-conditioned matrix
b7 = randn(100, 1);

% Time our implementation
tic;
x7_ours = mySolver(A7, b7);
time_ours = toc;

% Time MATLAB's implementation
tic;
x7_matlab = A7 \ b7;
time_matlab = toc;

% Compare results
diff7 = norm(x7_ours - x7_matlab);
residual7 = norm(A7*x7_ours - b7);

fprintf('  Our solver time: %.4f seconds\n', time_ours);
fprintf('  MATLAB solver time: %.4f seconds\n', time_matlab);
fprintf('  Speed ratio (ours/MATLAB): %.2f\n', time_ours/time_matlab);
fprintf('  Solution difference: %.2e\n', diff7);
fprintf('  Solution residual: %.2e\n', residual7);

fprintf('\n=== All tests completed ===\n');