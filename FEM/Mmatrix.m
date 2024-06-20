function is_m_matrix = Mmatrix(A)
    % Check if matrix A is a square matrix
    [m, n] = size(A);
    if m ~= n
        error('Matrix A must be square.');
    end
    
    % Check if all off-diagonal elements are non-positive
    for i = 1:n
        for j = 1:n
            if i ~= j && A(i, j) > 0
                is_m_matrix = false;
                return;
            end
        end
    end
    
    % Check if the matrix can be written as A = sI - B where B >= 0
    % and s >= rho(B)
    % Since B = sI - A, we need to find s
    % We can choose s as the maximum diagonal element of A
    s = max(diag(A));
    B = s * eye(n) - A;
    
    if any(B(:) < 0)
        is_m_matrix = false;
        return;
    end
    
    % Check if s >= rho(B)
    rho_B = max(abs(eig(B)));
    if s < rho_B
        is_m_matrix = false;
        return;
    end
    
    % Check if all eigenvalues of A have non-negative real parts
    eigenvalues_A = eig(A);
    if any(real(eigenvalues_A) < 0)
        is_m_matrix = false;
        return;
    end
    
    % If all checks pass, A is an M-matrix
    is_m_matrix = true;
end
