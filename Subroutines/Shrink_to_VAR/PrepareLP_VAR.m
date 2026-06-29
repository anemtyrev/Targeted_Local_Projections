function output = PrepareLP_VAR(input, H_max, H_min, P, one_matrix, which_irf_y, which_irf_x)
% PrepareLP_VAR: Prepare data for Local Projection estimation
%
% Handles three cases:
%   1. which_irf_y ~= which_irf_x: Response and shock are different variables
%   2. which_irf_y == which_irf_x: Own-response (e.g., G -> G)
%
% Inputs:
%   input       - struct with field y_t (T x k matrix of variables)
%   H_max       - maximum horizon
%   H_min       - minimum horizon (usually 0)
%   P           - number of lags
%   one_matrix  - 1 for stacked matrix, 0 for cell array by horizon
%   which_irf_y - index of response variable
%   which_irf_x - index of shock variable

if nargin <= 5
    which_irf_y = 1;
    which_irf_x = 2;
end

P = P - 1;
y_raw = input.y_t;
k = size(y_raw, 2);

% === Handle variable ordering based on whether own-response or not ===
if which_irf_y == which_irf_x
    % OWN-RESPONSE CASE: shock and response are the same variable
    which_y = y_raw(:, which_irf_y);
    which_x = which_y;  % Same variable
    
    % Other variables (excluding the shock/response variable)
    other = 1:k;
    idx_other = (other ~= which_irf_y);
    
    % Reorder: [response/shock, other variables]
    y = [which_y y_raw(:, idx_other)];
    
    % Create lagged controls
    w = lagmatrix(y, 1:P);
    newData = cat(2, y, w);
    newData(any(isnan(newData), 2), :) = [];
    T = size(newData, 1);
    
    % Extract variables
    y = newData(:, 1);              % Response variable (same as shock)
    x = newData(:, 1);              % Shock variable (same as response)
    w = newData(:, 2:size(newData, 2));  % All other contemporaneous + all lags
    
else
    % STANDARD CASE: shock and response are different variables
    which_y = y_raw(:, which_irf_y);
    which_x = y_raw(:, which_irf_x);
    
    % Other variables (excluding both shock and response)
    other = 1:k;
    idx_other = (other ~= which_irf_y & other ~= which_irf_x);
    
    % Reorder: [response, shock, other variables]
    y = [which_y which_x y_raw(:, idx_other)];
    
    % Create lagged controls
    w = lagmatrix(y, 1:P);
    newData = cat(2, y, w);
    newData(any(isnan(newData), 2), :) = [];
    T = size(newData, 1);
    
    % Extract variables
    y = newData(:, 1);              % Response variable
    x = newData(:, 2:k);            % Shock + other contemporaneous variables
    w = newData(:, k+1:size(newData, 2));  % All lags
end

% === Build LP matrices ===
if one_matrix == 1
    % STACKED MATRIX VERSION: creates 1 matrix with dimensions TH x PH
    w = [ones(size(w, 1), 1) w];
    H_max_new = H_max + 1 - H_min;
    
    idx = nan((H_max + 1) * T, 2);
    Y   = nan((H_max + 1) * T, 1);
    Xb  = zeros((H_max + 1) * T, H_max_new);
    Xc  = zeros((H_max + 1) * T, H_max_new, size(w, 2));
    
    for t = 1:T - H_min
        idx_beg = (t - 1) * H_max_new + 1;
        idx_end = t * H_max_new;
        idx(idx_beg:idx_end, 1) = t;
        idx(idx_beg:idx_end, 2) = H_min:H_max;
        
        % y
        y_range = (t + H_min):min((t + H_max), T);
        Y(idx_beg:idx_end) = [y(y_range); nan(H_max_new - length(y_range), 1)];
        
        % x - handle scalar vs vector
        if size(x, 2) == 1
            Xb(idx_beg:idx_end, :) = eye(H_max_new) * x(t);
        else
            Xb(idx_beg:idx_end, :) = eye(H_max_new) * x(t, 1);  % First column is shock
        end
        
        for i = 1:size(w, 2)
            Xc(idx_beg:idx_end, :, i) = eye(H_max_new) * w(t, i);
        end
    end
    
    X = Xb;
    Xmat = ones(length(Xc), 1);
    for i = 1:size(w, 2)
        Xmat = [Xmat Xc(:, :, i)];
    end
    Xmat(:, 1) = [];
    
    select = isfinite(Y);
    idx = idx(select, :);
    Y   = Y(select);
    X   = X(select, :);
    Xmat = Xmat(select, :);
    
    % Project out controls
    Xpr = X - Xmat * (Xmat \ X);
    Ypr = Y - Xmat * (Xmat \ Y);
    X = Xpr;
    Y = Ypr;
    
else
    % CELL ARRAY VERSION: each horizon in a separate cell
    X   = cell(H_max - H_min + 1, 1);
    Y   = cell(H_max - H_min + 1, 1);
    idx = cell(H_max - H_min + 1, 1);
    
    for h = H_min:H_max
        Y_horizon = lagmatrix(y, -h);        % Create Y vector h-horizons in the future
        select    = isfinite(Y_horizon);     % Selection vector for non-nan observations
        Y_horizon = Y_horizon(select);       % Select Y
        
        % Handle x depending on whether own-response or not
        if size(x, 2) == 1
            % Own-response case or single shock variable
            X_horizon = [x(select) w(select, :)];
        else
            % Standard case: x contains shock + other contemporaneous
            X_horizon = [x(select, :) w(select, :)];
        end
        
        % Populate the cells
        X{h - H_min + 1} = X_horizon;
        Y{h - H_min + 1} = Y_horizon;
        idx{h - H_min + 1} = [(1:length(Y_horizon))' h * ones(length(Y_horizon), 1)];
    end
end

% Store data
output.X = X;
output.Y = Y;
output.idx = idx;

end