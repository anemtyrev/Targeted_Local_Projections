% BootstrapConfInt_SLP
%
% This function performs a bootstrap procedure to estimate confidence intervals for
% the coefficients of a Smooth Local Projection (SLP). The function generates bootstrapped
% samples of the coefficients using a block bootstrap resampling technique and computes
% the confidence intervals using both the standard and studentized bootstrap methods.
%
% Syntax:
%    output = BootstrapConfInt_SLP(obj, lambda, alpha, bootstrapN, block_length)
%
% Inputs:
%    obj           - (struct) An object containing the following fields:
%                    - X (matrix): The independent variable matrix (size: T x k).
%                    - Z (matrix): The instrument matrix (if required by your model).
%                    - C (matrix): The matrix related to the penalty term.
%                    - Y (vector): The dependent variable vector (size: T x 1).
%                    - h (integer): The maximum horizon for forecasting.
%                    - H_min (integer): The minimum horizon for forecasting.
%                    - idx (matrix): Indices representing time periods and corresponding horizons.
%                    - nlag (integer): The number of lags in the model.
%    lambda        - (scalar) A penalty term for regularization (used in the coefficient estimation).
%    alpha         - (scalar) The significance level for constructing the confidence intervals (e.g., 0.05 for 95% confidence intervals).
%    bootstrapN    - (integer) The number of bootstrap samples to generate.
%    block_length  - (integer) The block length used for the bootstrap resampling.
%
% Outputs:
%    output        - (struct) A structure containing the following fields:
%                    - conf_int (matrix): A 2 x (h_max - h_min + 1) matrix containing the lower and upper bounds
%                      of the confidence intervals (from the standard bootstrap).
%                    - conf_int2 (matrix): A 2 x (h_max - h_min + 1) matrix containing the lower and upper bounds
%                      of the confidence intervals (from the studentized bootstrap).
%                    - beta_sort (matrix): The sorted bootstrap coefficient estimates across all iterations.
%
% Description:
% This function uses a block bootstrap procedure to resample the residuals and estimate the
% coefficients in a penalized linear regression model (SLP). The bootstrapped samples are then
% used to compute confidence intervals for the model coefficients, both using the standard
% bootstrap and studentized bootstrap approaches. The block size for the resampling is determined
% by `block_length`, and the number of bootstrap iterations is specified by `bootstrapN`.
%
%
% Notes:
% - The block bootstrap ensures that the residuals are resampled in blocks to preserve the time-series dependencies.
% - The studentized bootstrap is applied to adjust for heteroskedasticity in the errors by normalizing the coefficients with their standard errors.
% - The function provides two types of confidence intervals: the standard and the studentized.
% - The function performs a sanity check to ensure that the bootstrapped residuals (`U_boot`) match the expected size.


function output = BootstrapConfInt_SLP(obj,lambda,alpha,bootstrapN,use_studentized_boot,use_OLS_beta_for_resampling)

    if nargin <=5
        use_OLS_beta_for_resampling =1;
    end    

    %Unpacking the structure
    X        = obj.X;
    Z        = obj.Z;
    C        = obj.C;
    Y        = obj.Y;
    h_max    = obj.h;
    h_min    = obj.H_min;
    idx      = obj.idx;
    nlag     = obj.nlag;
    T        = max(idx(:,1));
    
    %running OLS
    if use_OLS_beta_for_resampling == 1
        beta = inv(X'*X)*X'*Y;
    else
        beta = inv(X'*X + lambda*C)*X'*Y;
    end
    U    = Y - X*beta;     
    block_length_T = ceil(T^(1/3));
    number_of_blocks = ceil(T / block_length_T);  % Calculate the number of blocks

    
    %creating empty array, first observation is observed beta
    beta_boot = nan(bootstrapN,h_max-h_min+1);  
    beta_boot(1,1:h_max-h_min+1) = beta(1:h_max-h_min+1);
    if use_studentized_boot == 1
        HAC_boot{1} = var_cov_quadratic(obj,lambda,alpha/2,1);
        VAR_boot(1,:) = sqrt(diag(HAC_boot{1}.VAR(1:h_max-h_min+1,1:h_max-h_min+1)));
    end
    % Bootstrap
    for b=2:bootstrapN   
        
        % Generate overlapping block start points
        start_points = randi(T - h_max - h_min-block_length_T, number_of_blocks, 1);  % Valid range for block starts
        block_indices = bsxfun(@plus, start_points, 0:block_length_T-1);    % Each row is a block
        indx_time_pool = reshape(block_indices', [], 1);  % Resampling pool of time indices
        indx_time_pool = indx_time_pool(1:T);             % Ensure we get exactly T samples
        
        % Initialize bootstrapped residual vector
        U_boot = [];
        
        % Loop over each t = 1 to T 
        for t = 1:T
            rows_t = find(idx(:,1) == t);         % Rows corresponding to original time t
            n_h_t = length(rows_t);               % Number of horizons at time t
        
            % Directly take a time index from indx_time_pool 
            t_star = indx_time_pool(t);   % Take the t-th index from the pool (overlapping blocks)
            rows_star = find(idx(:,1) == t_star);   % Rows at resampled time t_star
        
            % Ensure enough horizons are available in the resampled block
            if length(rows_star) >= n_h_t
                % Append the truncated residuals
                U_boot = [U_boot; U(rows_star(1:n_h_t))];
            else
                warning('Not enough horizons for time %d with resampled time %d', t, t_star);
            end
        end
        
        assert(size(U_boot, 1) == size(Y, 1), ...
            'Error: The size of U_boot does not match the number of rows in Y_boot!'); 
        
        % Subtract mean of U_boot per horizon using idx(:,2)
        U_boot_demeaned = U_boot;  % Initialize
        unique_horizons = unique(idx(:,2));

        for h = unique_horizons'
            rows_h = (idx(:,2) == h);                       % Logical index for horizon h
            mean_h = mean(U_boot(rows_h));                 % Mean for horizon h
            U_boot_demeaned(rows_h) = U_boot(rows_h) - mean_h;
        end

        % Generate the bootstrapped dependent variable
        Y_boot = X * beta + U_boot_demeaned;
        % Y_boot = X * beta + U_boot - mean(U_boot);
        
        beta_temp = inv(X'*X + lambda*C)*X'*Y_boot;
        beta_boot(b,:) = beta_temp(1:h_max+1-obj.H_min);

        obj2 = [];
        obj2 = obj;
        obj2.Y = Y_boot;
        
        if use_studentized_boot == 1
            HAC_boot{b} = var_cov_quadratic(obj2,lambda,alpha/2,1);
            VAR_boot(b,:) = sqrt(diag(HAC_boot{b}.VAR(1:h_max-h_min+1,1:h_max-h_min+1)));
        end
        % if mod(b,5)==0
        %     disp(['iteration ' num2str(b)]);
        % end
        % send(queue, b);
    end

    %sort the betas and confidence intervals
    beta_sort = sort(beta_boot);
    conf_int(1,:) = beta_sort(round(0.01*bootstrapN*(alpha)),:);
    conf_int(2,:) = beta_sort(round(0.01*bootstrapN*(100-alpha)),:);

    %store output
    output.conf_int = conf_int';
    % output.beta_sort = beta_sort;
    output.variance = var(beta_sort);

    %studentised bootstrap CI
    if use_studentized_boot == 1    
        student_beta = (beta_boot(1,:)-beta_boot)./(VAR_boot);
        t_sort_2 = sort(student_beta);
        conf_int2(1,:) = beta_boot(1,:) + VAR_boot(1,:).*t_sort_2(round(0.01*bootstrapN*(alpha)),:);
        conf_int2(2,:) = beta_boot(1,:) + VAR_boot(1,:).*t_sort_2(round(0.01*bootstrapN*(100-alpha)),:);
        
        output.conf_int = conf_int2';
    end
end


