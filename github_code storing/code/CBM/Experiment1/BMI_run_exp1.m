%BMI_RUN_EXP1  Run hierarchical Bayesian model comparison for Experiment 1.
%
%   This script performs a model comparison across a set of candidate
%   behavioural models using the CBM Hierarchical Bayesian Inference (HBI)
%   procedure.  It assumes that the Laplace approximations of each model
%   (lap_*.mat) have already been produced by BMI_fit_models.m.  The script
%   executes the following steps:
%
%     1) Load the eight‑item Likert ratings from the Excel file and format
%        them into a CBM data cell array.
%     2) Specify the candidate models and corresponding Laplace fit files.
%     3) Run cbm_hbi to compute group‑level model evidences and subjects'
%        responsibilities.
%     4) Apply a null correction (cbm_hbi_null) to obtain protected
%        exceedance probabilities.
%     5) Compute pooled and individual R² statistics to assess predictive
%        accuracy using get_mean_r_squared_bmi and
%        get_individual_mean_r_squared_bmi.

clear; clc;

% -------- Load Experiment 1 data --------
xlsx_file = '/Users/foxxis/Documents/social induction/experiment data/exp1/output100.xlsx';
raw_tbl   = readtable(xlsx_file);
vars      = {'d1','d2','d3','d4','d5','d6','d7','d8'};
tbl       = raw_tbl(:, vars);

Y    = table2array(tbl);
N    = size(Y, 1);
data = cell(N, 1);
for s = 1:N
    data{s} = struct('y', Y(s,:));
end

% -------- Define models and fitted file names --------
models    = {@mentalizing_model, @level0_model, @level1_model, @level2_model, @baseline_model};
fcbm_maps = {'lap_mentalizing.mat', 'lap_level0.mat', 'lap_level1.mat', 'lap_level2.mat', 'lap_baseline.mat'};

% Name of the HBI output file
hbi_file = 'hbi_compare.mat';

% (1) Run HBI (this generates/overwrites hbi_file)
cbm_hbi(data, models, fcbm_maps, hbi_file);

% (2) Apply null correction for protected exceedance probabilities
cbm_hbi_null(data, hbi_file);

% (3) Load the corrected HBI results
S   = load(hbi_file, 'cbm');
cbm = S.cbm;

% Extract posterior quantities (for inspection)
pxp  = cbm.output.protected_exceedance_prob;
xp   = cbm.output.exceedance_prob; %#ok<NASGU>
resp = cbm.output.responsibility;

% -------- Compute R² metrics --------

% Best model per subject (highest responsibility)
[~, bestModelPerSubj] = max(resp, [], 2);
[mean_r2_best, r2_best] = get_mean_r_squared_bmi(1000, data, cbm, bestModelPerSubj);
fprintf('Best‑model‑per‑subject pooled mean r² = %.3f\n', mean_r2_best);

% Winning model for all subjects (highest protected exceedance probability)
[~, winningModel] = max(pxp);
winningModelPerSubj = repmat(winningModel, numel(data), 1);
[mean_r2_win, r2_win] = get_mean_r_squared_bmi(1000, data, cbm, winningModelPerSubj);
fprintf('Winning‑model‑for‑all pooled mean r² = %.3f\n', mean_r2_win);

% Individual‑level mean R² using best model per subject
[~, bestModelPerSubj] = max(resp, [], 2);
[individual_mean_r2_best, individual_r2_indiv_boot_best] = get_individual_mean_r_squared_bmi(1000, data, cbm, bestModelPerSubj);
fprintf('Best‑model‑per‑subject (individual‑level) mean r² = %.3f\n', individual_mean_r2_best);

% Individual‑level mean R² using the winning model for all subjects
[~, winningModel] = max(pxp);
winningModelPerSubj = repmat(winningModel, numel(data), 1);
[individual_mean_r2_win, individual_r2_indiv_boot_win] = get_individual_mean_r_squared_bmi(1000, data, cbm, winningModelPerSubj);
fprintf('Winning‑model‑for‑all (individual‑level) mean r² = %.3f\n', individual_mean_r2_win);

% -------- Generate summary plots --------
% Determine best model per subject again (for plotting)
[~, bestModelPerSubj] = max(resp, [], 2);

% 1) Overall mean ratings vs predicted mean ratings
plot_overall_means(data, cbm, bestModelPerSubj);

% 2) Mean ratings by strategy group with predicted means
plot_choices_by_strategy(data, cbm, bestModelPerSubj);

% 3) Distributions of choices by strategy and trial with beta predictions
plot_choice_distribution_beta(data, cbm, bestModelPerSubj);