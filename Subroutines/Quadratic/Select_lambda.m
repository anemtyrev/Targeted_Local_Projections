function R = Select_lambda(lambda,Str,method)
    
    X        = Str.X;
    C        = Str.C;
    Y        = Str.Y;
    % T        = Str.T;
    % h        = Str.h;
    % idx      = Str.idx; 
    Sigma    = Str.Sigma;
    
    %Weighting matrix
    % W_hat = inv(Sigma);
    W_hat2 = (X'*X);

    %Shrinkage matrix
    V_of_lambda = inv(X'*X + lambda*C)*X'*X;

    %Shinked and unshrinked estimators
    Beta_OLS = inv(X'*X)*X'*Y;
    Beta_SR = V_of_lambda*Beta_OLS;
    
    % URE (modified and original)
    % I use W_hat as X'*X, in principle W_hat=eye(size(X'*X,1)) also works
    % R_var = (Beta_SR-Beta_OLS)'*W_hat*(Beta_SR-Beta_OLS)...
    %     + 2*trace(W_hat*V_of_lambda*Sigma);
    R_x_prime_x = (Beta_SR-Beta_OLS)'*W_hat2*(Beta_SR-Beta_OLS)...
        + 2*trace(W_hat2*V_of_lambda*Sigma);

    % % Cross validation
    S = X*inv(X'*X + lambda*C)*X'; 
    CV = sum( ( (Y - S*Y)./(1-diag(S)) ).^2 );
    
    if method == 1
        R = R_var;
    elseif method == 2
        R = R_x_prime_x;
    elseif method == 3 
        R = CV;
    else
        disp('unknown method')
    end    
end

