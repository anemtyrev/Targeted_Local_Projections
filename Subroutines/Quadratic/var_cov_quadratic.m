function output= var_cov_new(obj,lambda,alpha,kernel,undersmooth)
   
    if nargin<4, kernel = 1; end
    if nargin<5, undersmooth=0; end



    X        = obj.X;
    Z        = obj.Z;
    C        = obj.C;
    Y        = obj.Y;
    T        = obj.T;
    h        = obj.h;
    idx      = obj.idx;
    nlag     = obj.nlag;

    beta     = inv(X'*X + lambda*C)*X'*Y;
    U        = Y - X*beta;    
    npar     = length(beta);

    %Kernel weights (first is 0.5 bc variance is embedded in summation)
    if kernel ==1
        weights = [0.5, (1-(1:nlag)./(nlag+1)) ];  
    elseif kernel ==2
        lags = 0:nlag;
        x = lags/nlag;
        argQS = 6*pi*x/5;
        w1 = 3./(argQS.^2);
        w2 = (sin(argQS)./argQS)-cos(argQS);
        wQS = w1.*w2;
        wQS(x == 0) = 0.5;
        weights= wQS; 
    end

    V       = zeros( npar , npar );

    %% Bread
        W       = inv(Z'*Z)  ;
        XZ      = X'*Z*W*Z' ; 
        XXP     = XZ*X + lambda*C;
        if undersmooth==1
            XXP = XZ*X + 0.1*lambda*C;
        end
        theta   = XXP \ ( XZ*obj.Y );
        cow     = X'*Z;
        bread   = XXP^-1;
      
    %% MEAT
   
    for j = 0:nlag
        GplusGprime = zeros( npar , npar );
        for t = (j+1):(T-h-1) 
            S1 = Z( idx(:,1)==t , : )' * U( idx(:,1)==t );
            S2 = Z( idx(:,1)==(t-j) , : )' * U( idx(:,1)==(t-j) );
            GplusGprime = GplusGprime + S1 * S2' + S2 * S1';
        end
        V = V + weights(j+1) * GplusGprime;
    end
    meat = cow *W *V* W * cow' ;
    
    %% Sandwitch
    VAR = bread * meat * bread;

    %% Confidence interval
    conf_int = nan(npar,2);
    st_error = sqrt(diag(VAR(1:npar,1:npar)));                
    conf_int(:,1) = beta(1:npar) + st_error*norminv(0.01*alpha);
    conf_int(:,2) = beta(1:npar) + st_error*norminv(0.01*(100-alpha));

    %% Output

    output.VAR      = VAR;
    output.nlag     = nlag;
    output.conf_int = conf_int;
    output.U        = U;
end