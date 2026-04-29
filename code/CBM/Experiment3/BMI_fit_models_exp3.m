%BMI_FIT_MODELS_EXP3  Fit behavioural models to Experiment 3 data using CBM.
%
%   This script loads the Experiment 3 data, converts it to a CBM‑compatible
%   format, defines Laplace priors for each behavioural model, and fits
%   the models using cbm_lap.  The fitted results are saved into .mat
%   files.  The models themselves are implemented in separate function
%   files: mentalizing_model_exp3.m, level0_model_exp3.m,
%   level1_model_exp3.m, level2_model_exp3.m and baseline_model_exp3.m.
%
%   Adjust the paths to the CBM package and the data file as necessary.
%
%   See also BMI_RUN_EXP3.

clear; clc;

% -------- Configure paths --------
cbm_root = '/Users/foxxis/Documents/MATLAB/cbm-master';
addpath(genpath(cbm_root));
savepath;

% -------- Load Experiment 3 data --------
data_file = '/Users/foxxis/Documents/social induction/experiment data/exp3/final_data.xlsx';
raw_tbl   = readtable(data_file);
vars      = {'p1','p2','p3','p4','p5','p6','p7','p8','p9','p10','p11','p12','p13','p14','p15','p16'};
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
v = 1.8;  % prior variance suggested by the original script

% Fit mentalizing model
prior_ment = struct('mean', zeros(1,1), 'variance', v);
cbm_lap(data, @mentalizing_model_exp3, prior_ment, 'lap_mentalizing_exp3.mat');

% Fit level-0 model
prior_lv0 = struct('mean', zeros(2,1), 'variance', v);
cbm_lap(data, @level0_model_exp3, prior_lv0, 'lap_level0_exp3.mat');

% Fit level-1 (diversity heuristic) model
prior_lv1 = struct('mean', zeros(3,1), 'variance', v);
cbm_lap(data, @level1_model_exp3, prior_lv1, 'lap_level1_exp3.mat');

% Fit level-2 (colour) model
prior_lv2 = struct('mean', zeros(2,1), 'variance', v);
cbm_lap(data, @level2_model_exp3, prior_lv2, 'lap_level2_exp3.mat');

% Fit baseline model
prior_base = struct('mean', 0, 'variance', v);
cbm_lap(data, @baseline_model_exp3, prior_base, 'lap_baseline_exp3.mat');

% End of script