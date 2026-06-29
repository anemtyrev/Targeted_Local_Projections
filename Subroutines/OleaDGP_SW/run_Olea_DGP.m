
function output = run_Olea_DGP(str,p_DGP,scheme,specifcication_DGP,scale,garch_flag,bypass_Tscale)

if nargin <= 4
    scale = 1;
    garch_flag=0;
end
if nargin <= 6
    bypass_Tscale = 0;
end

% scheme = 'mpshock';  % mpshock, lshock, or mprecursive
addpath("OleaDGP_SW/")
sim_setscheme_sw
horzs = str.H_min:str.H_max;

%--------------------------------------------------------------------------
% Specifications
%--------------------------------------------------------------------------

% Note: "Specifications" contains the following fields:
%   - spec: worst (worst-case), estp (estimated p), fixp (fixed p)
%   - boot: Run bootstrap. 
%   - longT: true: Run T=2000. false: Run T=240.

exercise          = 'varma';
fields            = {'spec', 'boot', 'longT'};
Specifications    = cell2struct(cell(length(fields), 1), fields);
% Specifications(1) = struct('spec', 'worst', 'boot', false, 'longT', false);
% Specifications(2) = struct('spec', 'estp',  'boot', false, 'longT', false);
Specifications(1) = struct('spec', specifcication_DGP,  'boot', false, 'longT', false);
% Specifications(2) = struct('spec', 'fixp',  'boot', false, 'longT', true);

% Setup
load(['varma_sw_dgps_' scheme])  % Load dynare results (inputs/ folder is on the path via genpath)

spec  = Specifications(1).spec;
boot  = Specifications(1).boot;
longT = Specifications(1).longT;

sim_setup_general  % General setup
dgp.T = str.T+str.T_burn;
dgp.T_scale = 200;
sim_setup_sw
dgp.ps = p_DGP;% sw setup
dgps = dgp.ps;     % dgps considered by sw

numdgp  = size(dgps,2);                 % no. of DGPs
numhorz = length(settings.est.horzs);   % no. of horizons
numspec = length(specs);                % no. of regression specifications
numrep  = settings.sim.numrep;          % no. of repetitions

spec_shared = {'resp_ind',  settings.est.resp_ind, ...
               'innov_ind', settings.est.innov_ind, ...
               'alpha',     settings.est.alpha,...
               'no_const',  settings.est.no_const,...
               'se_homosk', settings.est.se_homosk,...
               'boot_num',  settings.est.boot_num};

dgp.irs_true       = nan(numdgp, numhorz);
dgp.var_asymp_covg = nan(numdgp, numhorz);

%% RUN SIMULATIONS

%----------------------------------------------------------------
% Placeholder Matrices
%----------------------------------------------------------------

estims = zeros(numdgp, numspec, numhorz, numrep);
ses    = estims;

cis_lower = zeros(numdgp, numspec, numhorz, 4, numrep);
cis_upper = cis_lower;

i_dgp = 1;
sim_run_setupdgp_sw

if bypass_Tscale
    dgp.zeta = 0;
end

dgp.alpha_tilde(2:end,:,:) = scale * dgp.T_scale^(dgp.zeta) * dgp.alpha_tilde(2:end,:,:); 
set_up_varma

data_y = generate_data(dgp,garch_flag);

output.y_t = data_y(str.T_burn+2:end,:);
output.dgp = dgp;

end