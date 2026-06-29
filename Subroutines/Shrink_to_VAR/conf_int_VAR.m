function output = conf_int_VAR(input,bootN,str,alpha,random_seed,index_boot,normal_boot) 

if nargin <=5
    index_boot  = 0;
    normal_boot = 1;
end   
%unpacking the structure
Y      = input.Y;         
X      = input.X;         
A_hat  = input.A_hat; 
% Y_boot = nan(size(Y,1),size(Y,2),bootN);           

% %create seed for the parfor loop
% random_seed = randi([0 (2^32-1)],bootN,1);

%dimentions of the data
[T,k] = size(Y);
p = length(A_hat)/k;
block_length = ceil(T^(1/3));
% block_length = 10;

%first observation
est_boot{1}      = input;
IRF(:,:,:,1)     = irf_function(input.A_hat_f,input.B_0_inv2,str);
% Y_boot(:,:,1)    = Y;

%create vector of residuals
%% 
U = Y - X*A_hat;
UU = input.resid(:,2);

%% check if LP resid are correct
        % aaa = PrepareLP_VAR(input,str.H_max,0,str.P_LP,0,2,1);
        % X2 = aaa.X{1}
        % Y2 = aaa.Y{1}
        % 
        % X1 = X2(:,1);
        % X_mat = [ones(size(X2,1),1) X2(:,2:end)];
        % 
        % Xpr = X1 - X_mat * (X_mat\X1);
        % Ypr = Y2 - X_mat * (X_mat\Y2);
        % 
        % Beta_OLS = Xpr\Ypr;
        % Eps = Ypr - Xpr*Beta_OLS';
        % 
        % X3 = Xpr; 
        % Y3 = Ypr;
        % 
        % resid = Y3 - X3 * (X3\Y3)
 %% bootstrap   

for b=2:bootN
    
    if index_boot==0
    %seed in the parfor loop
        % rng(random_seed(b-1))

        %resample data
        number_of_blocks = ceil(T/block_length);
        temp = randi(number_of_blocks-2,number_of_blocks,1)*block_length + 1;
        temp2 = temp + (0:block_length-1);
        indx_boot2 = reshape(temp2',block_length*number_of_blocks,1);
        indx_boot2(T+1:end,:) = [];
    else    
        indx_boot2 = index_boot(:,b);
    end    
    
    U_boot        = U(indx_boot2,:);
    U_boot_demean = U_boot - mean(U_boot);

    UU_boot       = UU(indx_boot2); 
    UU_boot       = UU_boot - mean(UU_boot);

    %recursive bootstrap 
    y_t = [];
    y_t(1:p,:) = U_boot_demean(1:p,:);
    for t=(p+1):T
        aaa = [];
        for i = 1:p
            aaa(i,:) = (A_hat(k*i - k+1:k*i,:)'*y_t(t-i,:)')';
        end
        y_t(t,:) = sum(aaa,1) + U_boot_demean(t,:);
    end    

    % residual bootstrap
    if normal_boot ==1
        Y_boot = y_t(:,:); 
        X_boot = lagmatrix(Y_boot,1:str.P_VARq);
    
        select = isfinite(X_boot);  
        Y_b = Y_boot(select(:,end),:);
        X_b = X_boot(select(:,end),:);
    
        %estimate on bootstrap data
        est_boot{b}       = estim_VAR(Y_b,str.P_VARq,X_b);
        IRF(:,:,:,b)      = irf_function(est_boot{b}.A_hat_f,est_boot{b}.B_0_inv2,str); 
        
        for i=str.H_min:str.H_max
            [Ci{i+1} Ci_bar{i+1}] = VAR_StructuralResponses(input.A_hat',input.Sigma_hat,i);
            variance2.tot(:,i+1,b) =  diag(Ci{i+1} * input.Sigma_alpha * Ci{i+1}' + Ci_bar{i+1}*input.Sigma_sigma*Ci_bar{i+1}' )/T;
            variance2.b0_part(:,i+1,b) =  diag(Ci_bar{i+1}*input.Sigma_sigma*Ci_bar{i+1}' )/T;
            variance2.A_part(:,i+1,b) =  diag(Ci{i+1} * input.Sigma_alpha * Ci{i+1}' )/T;
        end  


    else
        X2 = Y(:,1);
        Y2 = Y(:,2);

        yyy = Y(:,2);
        xxx = [Y(:,1) X ones(size(X2,1),1)];

        beta1 = xxx \ yyy;

        X1 = X2(:,1);
        X_mat = [ones(size(X2,1),1) X];

        Xpr = X1 - X_mat * (X_mat\X1);
        Ypr = Y2 - X_mat * (X_mat\Y2);

        Beta_OLS = Xpr\Ypr;

        Y_boot(1,1) = Y2(1,1);
        for i=str.P_LP:T
            Y_boot(i,1) = [X2(i,1) X(i,1) Y_boot(i-1,1) 1]*beta1 +  UU_boot(i,1);
        end    
        % Y_boot = Xpr*Beta_OLS' + UU_boot;
        Y2_boot = [Y(:,1) Y_boot];

        % mean(Y_boot)
        % std(Y_boot)
        % std(Y)

        X_boot = lagmatrix(Y2_boot,1:str.P_VARq);
    
        % beta = inv(X_b'*X_b)*X_b'*Y_b

        select = isfinite(X_boot);  
        Y_b = Y2_boot(select(:,end),:);
        X_b = X_boot(select(:,end),:);

        est_boot{b}       = estim_VAR(Y_b,str.P_VARq,X_b);
        IRF(:,:,:,b)      = irf_function(est_boot{b}.A_hat_f,est_boot{b}.B_0_inv2,str);

        if b==100
            sadada =232;
        end
        % A_0 = eye(k);
        % for i = 2:k
        %     Y1 = Y_boot(:,i);
        %     A  = Y(:,1:i-1); 
        %     X1 = [A X ones(T,1)];
        %     beta1 = X1\Y1;
        %     A_0(i,1:i-1) = -beta1(1:i-1)';
        %     resid(:,i) = Y1 - X1*beta1;
        %     % moment(:,i) = resid(:,i).*X1;
        % end
        % B_0 = inv(A_0);
        % IRF(:,:,:,b)      = irf_function(A_boot,B_0,str);
    end
end

%% sort and find quantile IRFs
irfs = reshape(permute(IRF,[2 1 3 4]),k^2,str.H_max-str.H_min+1,[]);
IRF_new  = sort(irfs,3);
CI_lower = IRF_new(:,:,round(0.01*bootN*(alpha/2)));
CI_upper = IRF_new(:,:,round(0.01*bootN*(100-alpha/2)));

% get delta method CI
for i=str.H_min:str.H_max
   [Ci{i+1} Ci_bar{i+1}] = VAR_StructuralResponses(input.A_hat',input.Sigma_hat,i);
    variance.tot(:,i+1) =  diag(Ci{i+1} * input.Sigma_alpha * Ci{i+1}' + Ci_bar{i+1}*input.Sigma_sigma*Ci_bar{i+1}' )/T;
    variance.b0_part(:,i+1) =  diag(Ci_bar{i+1}*input.Sigma_sigma*Ci_bar{i+1}' )/T;
    variance.A_part(:,i+1) =  diag(Ci{i+1} * input.Sigma_alpha * Ci{i+1}' )/T;
end  

%studentised bootstrap CI
student_beta = (irfs(:,:,1)-irfs)./(sqrt(variance2.tot));
t_sort_2 = sort(student_beta,3);
conf_int2(:,:,1) = irfs(:,:,1) + sqrt(variance.tot(:,:,1)).*t_sort_2(:,:,round(0.01*bootN*(alpha)));
conf_int2(:,:,2) = irfs(:,:,1) + sqrt(variance.tot(:,:,1)).*t_sort_2(:,:,round(0.01*bootN*(100-alpha)));

%delta_method CI
conf_int3(:,:,1) = irfs(:,:,1) - sqrt(variance.tot(:,:,1)).*norminv(0.01*(100-alpha/2));
conf_int3(:,:,2) = irfs(:,:,1) + sqrt(variance.tot(:,:,1)).*norminv(0.01*(100-alpha/2));

%store and reshape output
output.lower = CI_lower;
output.upper = CI_upper; 
output.upper_reshaped   = reshape(CI_upper,k^2,str.H_max-str.H_min+1);
output.lower_reshaped   = reshape(CI_lower,k^2,str.H_max-str.H_min+1);
output.conf_int         = [output.lower_reshaped(str.which_irf_var,:) ; output.upper_reshaped(str.which_irf_var,:)]';
output.conf_int_student = reshape(conf_int2(str.which_irf_var,:,:),str.H_max-str.H_min+1,[]);
output.conf_int_delta   = reshape(conf_int3(str.which_irf_var,:,:),str.H_max-str.H_min+1,[]);
output.var_delta        = variance.tot;
output.standard_dev_est = std(IRF_new,[],3);
output.standard_dev_est = std(IRF_new,[],3);
output.Y_boot           = Y_boot;
% output.est_boot         = est_boot;
output.median_beta      = median(IRF_new,3);
output.mean_beta        = mean(IRF_new,3);
% output.beta             = IRF(:,:,:,:);

end