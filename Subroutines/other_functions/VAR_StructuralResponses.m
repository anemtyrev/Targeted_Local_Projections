function [Ci, Ci_bar] = VAR_StructuralResponses(A, Sigma_u, i)
    % VAR_StructuralResponses calculates Gi and H matrices for structural VAR.
    %
    % Inputs:
    %   A         - VAR coefficient matrices [K x Kp]
    %   p         - Lag order of the VAR model
    %   Sigma_u   - Covariance matrix of residuals [K x K]
    %   i         - Horizon of interest
    %
    % Outputs:
    %   Gi        - Matrix G_i for structural impulse responses at horizon i
    %   H         - Matrix H for Cholesky decomposition adjustment

    % Dimensions
    K = size(Sigma_u, 1);  % Number of variables
    p = size(A,2)/K;
    
    % Step 1: Cholesky decomposition of Sigma_u
    B0_inv = chol(Sigma_u, 'lower'); 
    
    % Step 2: Define the selection matrix J
    J = [eye(K), zeros(K, K * (p - 1))]; % Selection matrix of size K x Kp

    % Step 3: Construct the companion matrix A_comp
    A_comp = [A; [eye(K * (p - 1))  zeros(K * (p - 1), K)]];
    
    % Step 4: Compute Gi for horizon i
    Gi = zeros(K^2,K^2*p);
    for m = 0:(i - 1)
        Am = A_comp'^(i - 1 - m);
        aaa = J *  Am;
        Psi_m = J * A_comp^m * J';
        Gi = Gi + kron(aaa, Psi_m);
    end

    % Step 5: Define elimination matrix L_K
    L_K = eliminationMatrix(K);
    
    % Step 6: Define commutation matrix K_KK
    K_KK = commutationMatrix(K, K);
    
    % Step 7: Compute H matrix
    H = L_K' * inv(L_K * (eye(K^2) + K_KK) * (kron(B0_inv, eye(K))) * L_K');

    Ci = kron(B0_inv', eye(K)) * Gi;

    Theta_i = J * A_comp^(i) * J'; % Theta_i = J * A^i * J'
    Ci_bar = kron(eye(K), Theta_i) * H;


end

% Helper function: Elimination matrix
function Lm = eliminationMatrix(m)
    % Creates the elimination matrix for symmetric matrices of size m x m
    Lm = zeros(m * (m + 1) / 2, m^2);
    idx = 1;
    for col = 1:m
        for row = col:m
            linearIdx = (row - 1) * m + col;
            Lm(idx, linearIdx) = 1;
            idx = idx + 1;
        end
    end
end

% Helper function: Commutation matrix
function Kmn = commutationMatrix(m, n)
    % Creates the commutation matrix of size mn x mn
    Kmn = zeros(m * n, m * n);
    for i = 1:m
        for j = 1:n
            Kmn((j - 1) * m + i, (i - 1) * n + j) = 1;
        end
    end
end
