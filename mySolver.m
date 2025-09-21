function x = mySolver(A, b)
%% mySolver - Linear System Solver using LU Decomposition
% This function solves the linear system Ax = b using LU decomposition
% with partial pivoting for improved numerical stability.
%
% Input:
%   A - Coefficient matrix (n x n)
%   b - Right-hand side vector (n x 1) or matrix (n x m)
%
% Output:
%   x - Solution vector (n x 1) or matrix (n x m)
%
% Method:
%   1. Decompose A = PLU using myLUDecomposer
%   2. Solve Ly = Pb using forward substitution
%   3. Solve Ux = y using back substitution
%
% Requirements:
%   - Matrix A must be square and non-singular
%   - Vector/matrix b must have compatible dimensions
%
% Author: Linear System Solver Implementation
% Date: Current

%% Input Validation
% Check if A is a matrix
if ~ismatrix(A) || ~isnumeric(A)
    error('mySolver:InvalidMatrixA', 'Input A must be a numeric matrix.');
end

% Check if A is square
[m, n] = size(A);
if m ~= n
    error('mySolver:NotSquare', 'Coefficient matrix A must be square.');
end

% Check if b is numeric
if ~isnumeric(b)
    error('mySolver:InvalidVectorB', 'Input b must be numeric.');
end

% Check dimensions compatibility
[b_rows, b_cols] = size(b);
if b_rows ~= n
    error('mySolver:IncompatibleDimensions', ...
          'Dimensions of A (%d x %d) and b (%d x %d) are incompatible.', ...
          m, n, b_rows, b_cols);
end

% Check for NaN or Inf values
if any(any(~isfinite(A))) || any(any(~isfinite(b)))
    error('mySolver:NonFiniteValues', 'Input matrices contain NaN or Inf values.');
end

%% Step 1: LU Decomposition with Partial Pivoting
try
    [P, L, U] = myLUDecomposer(A);
catch ME
    % Re-throw the error with additional context
    newME = MException('mySolver:DecompositionFailed', ...
                      'LU decomposition failed: %s', ME.message);
    newME = addCause(newME, ME);
    throw(newME);
end

%% Step 2: Forward Substitution - Solve Ly = Pb
% Compute Pb
Pb = P * b;

% Initialize solution for forward substitution
y = zeros(n, b_cols);

% Forward substitution: L is lower triangular with unit diagonal
for i = 1:n
    y(i, :) = Pb(i, :);
    for j = 1:i-1
        y(i, :) = y(i, :) - L(i, j) * y(j, :);
    end
    % Note: L(i,i) = 1, so no division needed
end

%% Step 3: Back Substitution - Solve Ux = y
% Initialize solution for back substitution
x = zeros(n, b_cols);

% Back substitution: U is upper triangular
for i = n:-1:1
    x(i, :) = y(i, :);
    for j = i+1:n
        x(i, :) = x(i, :) - U(i, j) * x(j, :);
    end
    
    % Check for division by zero (should not happen if LU decomposition succeeded)
    if abs(U(i, i)) < eps * norm(A, 'fro')
        error('mySolver:SingularMatrix', ...
              'Matrix is singular during back substitution at position (%d,%d).', i, i);
    end
    
    x(i, :) = x(i, :) / U(i, i);
end

%% Optional Verification (can be commented out for production)
% Check residual for numerical verification
residual = norm(A * x - b, 'fro') / norm(b, 'fro');
if residual > 1e-10
    warning('mySolver:LargeResidual', ...
            'Large residual detected (%.2e). Solution may be inaccurate.', residual);
end

end