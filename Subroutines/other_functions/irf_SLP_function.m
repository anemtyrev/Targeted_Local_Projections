function output = irf_SLP_function(data,str,method,calculate_LP,calculate_SLP,use_studentized_boot,use_OLS_beta_for_resampling)

    if nargin <=5
        use_studentized_boot = 0;
        use_OLS_beta_for_resampling = 1;
    end  

    %Quadratic penalty
    C = Create_C(str.H_max+1-str.H_min,1,0);

    % Reshape the data
    data_reshaped = PrepareLP_VAR(data, str.H_max, str.H_min, str.P_LP, 1, str.which_irf_y,str.which_irf_x);

    % Get the Y and X matrix
    X        = data_reshaped.X;
    Y        = data_reshaped.Y;
    idx      = data_reshaped.idx;
    Beta_OLS = X\Y;
        
    if calculate_SLP==1

        T=str.T;

        nlag = str.H_max-str.H_min +1;
    
        %Create structure 
        Struct          = [];
        Struct.C        = C;
        Struct.h        = str.H_max;
        Struct.H_min    = str.H_min;
        Struct.T        = T;
        Struct.X        = X;
        Struct.Y        = Y;
        Struct.Z        = X;           %no instruments yet
        Struct.idx      = idx;
        Struct.theta    = Beta_OLS;
        Struct.S        = inv(X'*X);
        Struct.nlag     = nlag;
        Struct.AR       = data;
    
        % HAC estimator
        HAC = var_cov_quadratic(Struct,0,str.alpha);
        Struct.Sigma = HAC.VAR;
    
        %optimization options
        opts = optimset('Display','off','TolCon',1E-10,'TolFun',1E-10,'TolX',1E-10);
        lambda_start = 500;
        x_L = 0;
    
        % choose between 1 = 'R_NEW', 2 = 'R_OLD', 3 = 'CV'
        objfun = @(lambda)Select_lambda(lambda,Struct,method); 
        lambda_opt = fmincon(objfun,lambda_start,[],[],[],[],x_L,[],[],opts);
    
        Beta_SR_opt = inv(X'*X + lambda_opt*C)*X'*Y;
    end
    %% CI

    if calculate_SLP==1
        %Bootstrap Quadratic
        conf_int_boot = BootstrapConfInt_SLP(Struct,lambda_opt,str.alpha/2,...
            str.BootN,use_studentized_boot,use_OLS_beta_for_resampling);
        %HAC SLP
        HAC_shrink = var_cov_quadratic(Struct,lambda_opt,str.alpha/2,str.HAC_kernel);
        HAC_undersmooth = var_cov_quadratic(Struct,lambda_opt,str.alpha/2,str.HAC_kernel,1);

        %output
        output.SLP.irf                  = Beta_SR_opt;
        output.SLP.CI_HAC               = HAC_shrink;
        output.SLP.CI_HAC_UNDERSMOOTH   = HAC_undersmooth;
        output.SLP.CI_Boot              = conf_int_boot;
        output.SLP.lambda               = lambda_opt;

    end

    if calculate_LP==1
        %Bootstrap OLS
        conf_int_OLS = BootstrapConfInt_SLP(Struct,0,str.alpha/2,...
            str.BootN,use_studentized_boot);
        %HAC OLS
        HAC_OLS    = var_cov_quadratic(Struct,0,str.alpha/2,str.HAC_kernel);
            
        %output
        output.LP.irf     = Beta_OLS;
        output.LP.CI_HAC  = HAC_OLS;
        output.LP.CI_Boot = conf_int_OLS;
    end

end