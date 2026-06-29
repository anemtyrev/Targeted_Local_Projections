function data_y = generate_data(dgp,garch_flag)

    if nargin <= 1
        garch_flag=0;
    end
    

    % Simulate data ABCD system
    
    % Inputs:
    % A         n_s x n_s       state evolution
    % B         n_s x n_e       state shock response
    % A         n_y x n_s       observables evolution
    % D         n_y x n_e       observables shock response
    % T         1 x 1           sample size
    
    % Output:
    % data_y    (T+1) x 1       simulated data (with initial condition)


    % preparations
    
    T = dgp.T;
    
    A = dgp.ABCD.A;
    B = dgp.ABCD.B;
    C = dgp.ABCD.C;
    D = dgp.ABCD.D;
    
    n_s   = dgp.n_s;
    n_eps = dgp.n_eps;
    n_y   = dgp.n_y;
    
    % draw shocks
    
    if garch_flag == 1
        omega = 0.05;
        alpha = 0.10;
        beta  = 0.85;
        
        data_eps = zeros(T, n_eps);
        sigma2   = zeros(T, n_eps);
        raw_eps  = randn(T, n_eps);
        
        sigma2(1,:) = omega / (1 - alpha - beta);
        data_eps(1,:) = sqrt(sigma2(1,:)) .* raw_eps(1,:);

        % create_data
        for t = 2:T
            sigma2(t,:) = omega + alpha * data_eps(t-1,:).^2 + beta * sigma2(t-1,:);
            data_eps(t,:) = sqrt(sigma2(t,:)) .* raw_eps(t,:);
        end
    else
        data_eps = randn(T,n_eps);
    end    
    
    % simulate states
    
    s = zeros(n_s,1);
    data_s = NaN(T,n_s);
    for t = 1:T
        s = A * s + B * data_eps(t,:)';
        data_s(t,:) = s';
    end
    
    % simulate observables

    data_y = NaN(T,n_y);
    data_y(1,:) = (D * data_eps(1,:)')';
    for t = 2:T
        data_y(t,:) = (C * data_s(t-1,:)' + D * data_eps(t,:)')';
    end

    data_y = [zeros(1,n_y);data_y];

end