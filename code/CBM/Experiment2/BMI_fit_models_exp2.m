%BMI_FIT_MODELS_EXP2  Fit behavioural models to Experiment 2 data using CBM.
%
%   This script loads the Experiment 2 data, converts it to a CBM‑compatible
%   format, defines Laplace priors for each behavioural model, and fits the
%   models using cbm_lap.  The fitted results are saved into .mat files.
%   The models themselves are implemented in separate function files:
%     mentalizing_model_exp2.m, level0_model_exp2.m, level1_model_exp2.m,
%     level2_model_exp2.m and baseline_model.m.
%
%   Adjust the path to the CBM package and the data file as necessary.
%
%   See also BMI_RUN_EXP2.

clear; clc;

% -------- Configure paths --------
cbm_root = '/Users/foxxis/Documents/MATLAB/cbm-master';
addpath(genpath(cbm_root));
savepath;

% -------- Load Experiment 2 data --------
data_file = '/Users/foxxis/Documents/social induction/experiment data/exp2/final_data.csv';
raw_tbl   = readtable(data_file);
vars      = {'p1','p2','p3','p4','p5','p6','p7','p8'};
missing = setdiff(vars, raw_tbl.Properties.VariableNames);
assert(isempty(missing), 'Missing columns in %s: %s', data_file, strjoin(missing, ', '));
tbl = raw_tbl(:, vars);
Y   = table2array(tbl);
N   = size(Y,1);

% Convert to CBM data format: a cell array of structs with field y
data = cell(N,1);
for s = 1:N
    data{s} = struct('y', Y(s,:));
end

% -------- Set priors --------
v = 1;  % prior variance

% Fit mentalizing model
prior_ment = struct('mean', zeros(1,1), 'variance', v);
cbm_lap(data, @mentalizing_model_exp2, prior_ment, 'lap_mentalizing_exp2.mat');

% Fit level-0 model
prior_lv0 = struct('mean', zeros(2,1), 'variance', v);
cbm_lap(data, @level0_model_exp2, prior_lv0, 'lap_level0_exp2.mat');

% Fit level-1 (diversity heuristic) model
prior_lv1 = struct('mean', zeros(3,1), 'variance', v);
cbm_lap(data, @level1_model_exp2, prior_lv1, 'lap_level1_exp2.mat');

% Fit level-2 (colour) model
prior_lv2 = struct('mean', zeros(2,1), 'variance', v);
cbm_lap(data, @level2_model_exp2, prior_lv2, 'lap_level2_exp2.mat');

% Fit baseline model
prior_base = struct('mean', 0, 'variance', v);
cbm_lap(data, @baseline_model, prior_base, 'lap_baseline_exp2.mat');

% End of script
