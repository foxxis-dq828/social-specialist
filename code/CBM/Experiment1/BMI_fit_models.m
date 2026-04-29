%BMI_FIT_MODELS  Fit behavioural models to Experiment 1 data using CBM.
%
%   This script fits several behavioural models to the eight‑item Likert data
%   collected in Experiment 1.  It uses the Hierarchical Bayesian Inference
%   framework provided by the CBM package (Piray & Daw) to perform a Laplace
%   approximation (cbm_lap) for each model.  The models themselves are
%   implemented in separate function files located in the same directory:
%   mentalizing_model.m, level0_model.m, level1_model.m, level2_model.m and
%   baseline_model.m.
%
%   The script performs the following steps:
%     1) Add the CBM package to the MATLAB path (edit cbm_root if needed).
%     2) Load the Excel file containing eight columns (d1..d8) of Likert
%        ratings and convert it into a cell array of subject structures
%        required by CBM (one struct per subject with field y).
%     3) Define a prior variance (v) for all model parameters.
%     4) Call cbm_lap for each model, saving the fitted results into
%        individual .mat files (e.g. lap_mentalizing.mat, lap_level0.mat).

% -------- Configure paths --------
cbm_root = '/Users/foxxis/Documents/MATLAB/cbm-master';
addpath(genpath(cbm_root));
savepath;  % optional: persist the path across sessions

% Verify the CBM functions are available
which('cbm_lap', '-all')

% -------- Load Experiment 1 data --------
xlsx_file = '/Users/foxxis/Documents/social induction/experiment data/exp1/output100.xlsx';
raw_tbl   = readtable(xlsx_file);
vars      = {'d1','d2','d3','d4','d5','d6','d7','d8'};

% ensure the required columns exist
missing = setdiff(vars, raw_tbl.Properties.VariableNames);
assert(isempty(missing), 'Missing columns in %s: %s', xlsx_file, strjoin(missing, ', '));

tbl = raw_tbl(:, vars);
Y   = table2array(tbl);  % N × 8 numeric array
N   = size(Y, 1);

% Convert to CBM data format: a cell array of structs with field y
data = cell(N, 1);
for s = 1:N
    data{s} = struct('y', Y(s,:));
end

% -------- Set priors --------
% The variance v controls the scale of the zero‑mean Laplace prior on each
% parameter.  A value of 1.8 is suggested by the CBM manual.
v = 1.8;

% Mentalizing model: one free parameter (log_ce)
prior_ment = struct('mean', zeros(1,1), 'variance', v);
cbm_lap(data, @mentalizing_model, prior_ment, 'lap_mentalizing.mat');

% Level‑0 model: two free parameters (log_ce, log_n)
prior_lv0 = struct('mean', zeros(2,1), 'variance', v);
cbm_lap(data, @level0_model, prior_lv0, 'lap_level0.mat');

% Level‑1 (diversity heuristic) model: three free parameters (log_ce, log_n, log_e)
prior_lv1 = struct('mean', zeros(3,1), 'variance', v);
cbm_lap(data, @level1_model, prior_lv1, 'lap_level1.mat');

% Level‑2/colour model: two free parameters (log_ce, log_n)
prior_lv2 = struct('mean', zeros(2,1), 'variance', v);
cbm_lap(data, @level2_model, prior_lv2, 'lap_level2.mat');

% Baseline model: zero free parameters (uniform over responses)
prior_base = struct('mean', 0, 'variance', v);
cbm_lap(data, @baseline_model, prior_base, 'lap_baseline.mat');

% End of script