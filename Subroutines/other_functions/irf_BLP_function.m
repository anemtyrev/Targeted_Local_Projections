function output = irf_BLP_function(data,str)

k = str.k;

hyperPars.isrw                 = false; %rw prior

%initialize BLP prior on presample
modelSpec.presample            = true;

%set BLP prior type
modelSpec.priorType            ='VAR'; %BLP only, alternatives: 'VAR', 'RW', 'DSGE'; 

%set VAR specification, lags, IRFs horizon
modelSpec.nVARlags             = str.P_VAR; %number of lags in VAR (& prior)
modelSpec.nBLPlags             = str.P_VAR; %number of lags in BLP
modelSpec.nLPlags              = str.P_VAR; %number of lags in LP

modelSpec.nHorizons            =str.H_max;
modelSpec.bandsCoverage        =100-str.alpha;

%choose identification scheme
modelSpec.identification       ='CHOL'; %alternatives: 'PSVAR', 'CHOL';

%declare shock variable & shock size
shockVar                       ='FFR';       
shockSize                      =1;           %percentage points


%-HYPERPRIORS SETTINGS----------------------------------------------------%

%hyperpriors initial values
%
hyperPars.lambda                 =.4;    %tightness of VAR coeffs prior (the higher lambda the closer to OLS)
hyperPars.lambdaC                =1e5;   %intercept (you want this to be large)
hyperPars.lambdaP                =.4;    %tightness of projCoeffs prior (same role as above)
%
hyperPars.miu                    =1;     %sum of coefficients prior (constraint multiplier)
hyperPars.theta                  =2;     %cointegration prior (constraint multiplier) 
hyperPars.alpha                  =2;     %lag decaying coeff for NIW prior


%set hyperpriors options (matches GLP fields); if you want default values
%(when available) set to empty [];
hyperPriorsOptions.hyperpriors   =true;                  %find optimal hyperparameters: NO default option
hyperPriorsOptions.Vc            =1e5;                   %variance of the VAR constant (default=1e6)
hyperPriorsOptions.MNalpha       =[];                    %lag decaying coeff of NIW prior (default=2)
hyperPriorsOptions.pos           =find(~hyperPars.isrw); %position of stationary variables
hyperPriorsOptions.MNpsi         =false;                 %residual variance univariate AR(1) std (default=hyperprior)
hyperPriorsOptions.noc           =false;                 %sum of coefficients prior: NO default option
hyperPriorsOptions.sur           =false;                 %cointegration prior: NO default option
hyperPriorsOptions.Fcast         =false;                 %build forecasts: NO default option
hyperPriorsOptions.hz            =modelSpec.nHorizons;   %max forecast horizon: NO default option
hyperPriorsOptions.mcmc          =false;                 %run metropolis-hasting algorithm: NO default option
hyperPriorsOptions.Ndraws        =1200;                  %default=20k
hyperPriorsOptions.Ndrawsdiscard =200;                   %default=10k
hyperPriorsOptions.MCMCconst     =1;                     %default=1
hyperPriorsOptions.MCMCfcast     =false;                 %store forecast at each MCMC draw (default=true)
hyperPriorsOptions.MCMCstorecoeff=false;                 %store coefficients at each MCMC draw (default=true)
hyperPriorsOptions.initialValues =hyperPars;             %see above

hyperPriorsOptions.priorType = modelSpec.priorType;

%sampling from paramenters distribution
GibbsOptions.iterations          = 1200;
GibbsOptions.burnin              = 200;
GibbsOptions.jump                = 1;

hyperPriorsOptions.GibbsOptions  =GibbsOptions;
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-LOAD DATA IN MODEL STRUCTURE--------------------------------------------%

%relevant sample
presample = max(30,round(str.T/6));


%use relevant observations
macroVarData                =data.y_t(presample+1:end,:);
macroVarDates               =1:presample;

macroVarPreData             =data.y_t(1:presample,:);
macroVarPreDates            =presample+1:str.T;


macroVarNames               =string(1:k);
macroVarLabels              =string(1:k);


%load data structure
dataStructure.data          =macroVarData;
dataStructure.dates         =macroVarDates;
dataStructure.preSdata      =macroVarPreData;
dataStructure.preSdates     =macroVarPreDates;

dataStructure.varname       =macroVarLabels;
dataStructure.varLongName   =macroVarNames;


%load shock variable and normalization in model structure 
modelSpec.shockSize                       = zeros(1,k)';
modelSpec.shockSize(str.which_irf_x)      = 1;
modelSpec.shockVar                        = modelSpec.shockSize;


%load data in model specification
modelSpec.dataStructure     =dataStructure;

%-------------------------------------------------------------------------%
%                   ESTIMATE IMPULSE RESPONSE FUNCTIONS                   %
%-------------------------------------------------------------------------%


% % %estimate LP & compute IRFs
% IRF_LP         =IRFlocalProj(modelSpec);
% % 
% % %estimate VAR & compute IRFs
% IRF_BVAR       =IRFbayesianNIW(modelSpec,hyperPriorsOptions);

%estimate BLP & compute IRFs
[~,IRF_BLP]        = evalc('IRFbayesianLocalProj(modelSpec,hyperPriorsOptions)');

output.irf     = IRF_BLP.irfs(:,str.which_irf_y);
output.CI      = [IRF_BLP.irfs_l(:,str.which_irf_y) IRF_BLP.irfs_u(:,str.which_irf_y)];
output.var     = IRF_BLP.st_dev(:,str.which_irf_y)'.^2;



