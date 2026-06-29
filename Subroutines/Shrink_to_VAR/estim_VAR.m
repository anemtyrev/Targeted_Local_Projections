function output = estim_VAR(str,p,X)

if X==0
    Y = str.y_t;
    %create lagged matrix
    Y_lag = lagmatrix(Y,1:p);
    
    %deleting nan observatoions
    Y_lag(1:p,:) = [];
    Y(1:p,:) = [];
else
    Y = str;
    %create lagged matrix
    Y_lag = X;
end

%figuring out correct dimentions
[T,k] = size(Y);

%for checks
Beta = inv(Y_lag'*Y_lag)*Y_lag'*Y;

%vectorising the data and finding A_hat
Y_vec = reshape(Y,k*T,1);
X_vec = kron(eye(k),Y_lag);
A_vec = X_vec\Y_vec;
A_hat = reshape(A_vec,k*p,k);

%Finding variance
U = Y - Y_lag*A_hat;
% U = U ./ std(U);
Sigma_hat = cov(U);
B_0_inv = chol(Sigma_hat,'lower');

Z = constructLaggedMatrix(Y, p);

Z'*Z;

Sigma_alpha_hat2 = kron(Sigma_hat,inv(Y_lag'*Y_lag/T));
Sigma_alpha_hat = kron(inv(Y_lag'*Y_lag),T*Sigma_hat);



A_0 = eye(k);
for i = 2:k;
    Y1 = Y(:,i);
    A  = Y(:,1:i-1); 
    X1 = [A Y_lag ones(T,1)];
    beta1 = X1\Y1;
    A_0(i,1:i-1) = -beta1(1:i-1)';
    resid(:,i) = Y1 - X1*beta1;
end

% aaa= B_0_inv
% aaa(1,1) = aaa(1,1)./B_0_inv(1,1);
% bbb=aaa;
% bbb(2,:) = aaa(2,:)./aaa(2,2);

[D, L, K] = select_matrix(k);
D_plus = inv(D'*D)*D';
Sigma_sigma_hat = 2*D_plus*kron(Sigma_hat,Sigma_hat)*D_plus';

%FAKE VAR(1)
A_hat_f   = [A_hat'; eye(k*(p-1)) zeros(k*(p-1),k)];

forecast_resid = Y - Y_lag * Beta;

%generating output  
output.Y             = Y;
output.X             = Y_lag;
output.A_hat         = A_hat; 
output.A_hat_f       = A_hat_f;
output.B_0_inv       = B_0_inv;
output.Sigma_hat     = Sigma_hat;
output.Sigma_alpha   = Sigma_alpha_hat;
output.Sigma_sigma   = Sigma_sigma_hat;
output.B_0_inv2      = inv(A_0);
output.resid         = resid;
output.y_t           = Y;
output.project_resid = forecast_resid;

end