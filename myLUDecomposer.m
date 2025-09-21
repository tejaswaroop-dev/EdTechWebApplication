function [P, L, U] = myLUDecomposer(A)
%% myLUDecomposer - LU Decomposition with Partial Pivoting
% This function performs LU decomposition of a square matrix A using
% partial pivoting to ensure numerical stability.
%
% Input:
%   A - Square matrix to be decomposed (n x n)
%
% Output:
%   P - Permutation matrix (n x n)
%   L - Lower triangular matrix with unit diagonal (n x n)
%   U - Upper triangular matrix (n x n)
%
% The decomposition satisfies: P*A = L*U
%
% Requirements:
%   - Matrix A must be square and non-singular
%   - Function implements partial pivoting for numerical stability
%
% Author: LU Decomposition Implementation
% Date: Current

%% Input Validation
% Check if A is a matrix
if ~ismatrix(A) || ~isnumeric(A)
    error('myLUDecomposer:InvalidInput', 'Input A must be a numeric matrix.');
end

% Check if A is square
[m, n] = size(A);
if m ~= n
    error('myLUDecomposer:NotSquare', 'Input matrix A must be square.');
end

% Check if A is empty
if isempty(A)
    error('myLUDecomposer:EmptyMatrix', 'Input matrix A cannot be empty.');
end

% Check for NaN or Inf values
if any(any(~isfinite(A)))
    error('myLUDecomposer:NonFiniteValues', 'Input matrix A contains NaN or Inf values.');
end

%% Initialize matrices
U = A;                          % Upper triangular matrix (will be modified)
L = eye(n);                     % Lower triangular matrix (unit diagonal)
P = eye(n);                     % Permutation matrix

%% LU Decomposition with Partial Pivoting
for k = 1:n-1
    % Find the pivot (largest absolute value in column k from row k to n)
    [~, pivot_row] = max(abs(U(k:n, k)));
    pivot_row = pivot_row + k - 1;  % Adjust index to full matrix
    
    % Check for singular matrix (pivot too small)
    if abs(U(pivot_row, k)) < eps * norm(A, 'fro')
        error('myLUDecomposer:SingularMatrix', ...
              'Matrix is singular or nearly singular (rank deficient).');
    end
    
    % Swap rows if necessary (partial pivoting)
    if pivot_row ~= k
        % Swap rows in U
        temp = U(k, :);
        U(k, :) = U(pivot_row, :);
        U(pivot_row, :) = temp;
        
        % Swap rows in P
        temp = P(k, :);
        P(k, :) = P(pivot_row, :);
        P(pivot_row, :) = temp;
        
        % Swap corresponding elements in L (only the computed part)
        if k > 1
            temp = L(k, 1:k-1);
            L(k, 1:k-1) = L(pivot_row, 1:k-1);
            L(pivot_row, 1:k-1) = temp;
        end
    end
    
    % Eliminate entries below the diagonal in column k
    for i = k+1:n
        % Calculate multiplier
        L(i, k) = U(i, k) / U(k, k);
        
        % Update row i of U
        U(i, :) = U(i, :) - L(i, k) * U(k, :);
    end
end

% Final check for singularity
if abs(U(n, n)) < eps * norm(A, 'fro')
    error('myLUDecomposer:SingularMatrix', ...
          'Matrix is singular or nearly singular (rank deficient).');
end

%% Verification (optional, can be commented out for production)
% Verify that P*A = L*U (within numerical tolerance)
if nargout == 0 || max(max(abs(P*A - L*U))) > 1e-12 * norm(A, 'fro')
    warning('myLUDecomposer:NumericalError', ...
            'LU decomposition may have numerical errors. Check condition number of A.');
end

end