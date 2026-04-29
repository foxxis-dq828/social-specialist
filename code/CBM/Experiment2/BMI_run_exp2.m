%BMI_RUN_EXP2  Run hierarchical Bayesian model comparison for Experiment 2.
%
%   This script performs HBI model comparison across the candidate models
%   fitted to Experiment 2 data.  It assumes that the Laplace fits have
%   been produced by BMI_FIT_MODELS_EXP2.m.  The script executes:
%     1) Load the Experiment 2 data.
%     2) Specify the models and their corresponding fit files.
%     3) Run cbm_hbi and apply the null correction.
%     4) Compute pooled and individual-level mean R² statistics.
%     5) Generate summary plots comparing actual ratings and model predictions.
%
%   See also BMI_FIT_MODELS_EXP2.

clear; clc;

% -------- Load Experiment 2 data --------
data_file = '/Users/foxxis/Documents/social induction/experiment data/exp2/final_data.csv';
raw_tbl   = readtable(data_file);
vars      = {'p1','p2','p3','p4','p5','p6','p7','p8'};
tbl       = raw_tbl(:, vars);
Y    = table2array(tbl);
N    = size(Y, 1);
data = cell(N,1);
for s = 1:N
    data{s} = struct('y', Y(s,:));
end

% -------- Define models and fitted file names --------
models    = {@mentalizing_model_exp2, @level0_model_exp2, @level1_model_exp2, @level2_model_exp2, @baseline_model};
fcbm_maps = {'lap_mentalizing_exp2.mat', 'lap_level0_exp2.mat', 'lap_level1_exp2.mat', 'lap_level2_exp2.mat', 'lap_baseline_exp2.mat'};

% Name of the HBI output file
hbi_file = 'hbi_compare_exp2.mat';

% (1) Run HBI
cbm_hbi(data, models, fcbm_maps, hbi_file);

% (2) Apply null correction
cbm_hbi_null(data, hbi_file);

% (3) Load the corrected HBI results
S   = load(hbi_file, 'cbm');
cbm = S.cbm;

% Extract posterior quantities
pxp  = cbm.output.protected_exceedance_prob;
resp = cbm.output.responsibility;

% -------- Compute R² metrics --------
[~, bestModelPerSubj] = max(resp, [], 2);
[mean_r2_best, r2_best] = get_mean_r_squared_bmi_exp2(1000, data, cbm, bestModelPerSubj);
fprintf('Best‑model‑per‑subject pooled mean r² = %.3f\n', mean_r2_best);

[~, winningModel] = max(pxp);
winningModelPerSubj = repmat(winningModel, numel(data), 1);
[mean_r2_win, r2_win] = get_mean_r_squared_bmi_exp2(1000, data, cbm, winningModelPerSubj);
fprintf('Winning‑model‑for‑all pooled mean r² = %.3f\n', mean_r2_win);

[~, bestModelPerSubj] = max(resp, [], 2);
[individual_mean_r2_best, individual_r2_indiv_boot_best] = get_individual_mean_r_squared_bmi_exp2(1000, data, cbm, bestModelPerSubj);
fprintf('Best‑model‑per‑subject (individual‑level) mean r² = %.3f\n', individual_mean_r2_best);

[~, winningModel] = max(pxp);
winningModelPerSubj = repmat(winningModel, numel(data), 1);
[individual_mean_r2_win, individual_r2_indiv_boot_win] = get_individual_mean_r_squared_bmi_exp2(1000, data, cbm, winningModelPerSubj);
fprintf('Winning‑model‑for‑all (individual‑level) mean r² = %.3f\n', individual_mean_r2_win);

% -------- Generate summary plots --------
% Determine best model per subject again (for plotting)
[~, bestModelPerSubj] = max(resp, [], 2);

% 1) Overall mean ratings vs predicted mean ratings
plot_overall_means_exp2(data, cbm, bestModelPerSubj);

[~, bestModelPerSubj] = max(resp, [], 2);
plot_choices_by_strategy_exp2(data, cbm, bestModelPerSubj);
plot_choice_distribution_beta_exp2(data, cbm, bestModelPerSubj);


