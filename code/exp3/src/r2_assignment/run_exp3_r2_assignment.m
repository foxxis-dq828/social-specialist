function run_exp3_r2_assignment(projectDir, data, cbm, n_simulations)
% Each iteration draws one categorical response for every participant and trial
% from the model-predicted choice distribution.
%
% Two R2 summaries are reported:
%   1. Pooled R2: concatenate all valid participant-by-trial observations and
%      simulated choices, then calculate corr(Yobs, Ysim)^2 once.
%   2. Individualized R2: calculate corr(Yobs, Ysim)^2 separately within each
%      participant, then average the participant-level R2 values.

addpath(projectDir);

addpath(fullfile(projectDir, 'src'));

% set n_simulations = 10000 as default
if nargin < 4 || isempty(n_simulations)
    n_simulations = 10000;
end

% calculate the model best predict each participant
responsibility = cbm.output.responsibility;
[~, best_model_per_subj] = max(responsibility, [], 2);
best_model_per_subj = best_model_per_subj(:);

% calculate the group best model
group_freq = cbm.output.model_frequency;
[~, group_winner] = max(group_freq);
group_winner_per_subj = repmat(group_winner, numel(best_model_per_subj), 1);

fprintf('Individualized assignment uses each participant''s maximum-responsibility model.\n');
fprintf('Group-winning assignment uses model %d for every participant.\n', group_winner);

% set the dataframe for calculating r2
n_participants = numel(data);
n_trials = numel(data{1}.y);
Yobs = nan(n_participants, n_trials);
for s = 1:n_participants
    Yobs(s, :) = data{s}.y(:)';
end

pooled_r2_individualized = nan(n_simulations, 1);
pooled_r2_group_winner = nan(n_simulations, 1);
individualized_r2_individualized = nan(n_simulations, 1);
individualized_r2_group_winner = nan(n_simulations, 1);

valid_mask = ~isnan(Yobs) & (Yobs >= 1) & (Yobs <= 6);

% Running two simulations using the same random seed
seeds = randi(2^31 - 1, n_simulations, 1);
for sim = 1:n_simulations
    if mod(sim, 1000) == 0
        fprintf('Categorical simulation iteration %d/%d\n', sim, n_simulations);
    end
    % simulate dataset using best model for each subject
    rng(seeds(sim));
    [~, Ysim_individualized] = simulate_dataset(cbm, best_model_per_subj, data);
    % simulate dataset using group winner model for all subjects
    rng(seeds(sim));
    [~, Ysim_group_winner] = simulate_dataset(cbm, group_winner_per_subj, data);

    % for each simulations, calculate individualized/pooled r2
    % calculate r2 by pooling all trials across all participants
    pooled_r2_individualized(sim) = pooled_r2(Yobs, Ysim_individualized, valid_mask);
    pooled_r2_group_winner(sim) = pooled_r2(Yobs, Ysim_group_winner, valid_mask);
    % calculate r2 separately for each participant, then average across participants.
    individualized_r2_individualized(sim) = mean(subject_r2(Yobs, Ysim_individualized, valid_mask), 'omitnan');
    individualized_r2_group_winner(sim) = mean(subject_r2(Yobs, Ysim_group_winner, valid_mask), 'omitnan');
end

pooled_r2_diff = pooled_r2_individualized - pooled_r2_group_winner;
individualized_r2_diff = individualized_r2_individualized - individualized_r2_group_winner;

valid_pooled_diff = pooled_r2_diff(~isnan(pooled_r2_diff));
valid_individualized_diff = individualized_r2_diff(~isnan(individualized_r2_diff));

pooled_diff_ci = prctile(valid_pooled_diff, [2.5 97.5]);
individualized_diff_ci = prctile(valid_individualized_diff, [2.5 97.5]);

r2_results = table( ...
    mean(pooled_r2_individualized, 'omitnan'), ...
    mean(pooled_r2_group_winner, 'omitnan'), ...
    mean(valid_pooled_diff, 'omitnan'), ...
    pooled_diff_ci(1), ...
    pooled_diff_ci(2), ...
    mean(individualized_r2_individualized, 'omitnan'), ...
    mean(individualized_r2_group_winner, 'omitnan'), ...
    mean(valid_individualized_diff, 'omitnan'), ...
    individualized_diff_ci(1), ...
    individualized_diff_ci(2), ...
    'VariableNames', { ...
        'pooled_r2_individualized', ...
        'pooled_r2_group_winner', ...
        'pooled_r2_diff', ...
        'pooled_ci_low', ...
        'pooled_ci_high', ...
        'mean_subject_r2_individualized', ...
        'mean_subject_r2_group_winner', ...
        'mean_subject_r2_diff', ...
        'mean_subject_ci_low', ...
        'mean_subject_ci_high', ...
    });

disp(r2_results);
end

% Helper functions
function [sim_data, Ysim] = simulate_dataset(cbm, modelChoicePerSubj, templateData)
%SIMULATE_DATASET Generate a synthetic CBM dataset from HBI parameters by sampling
% one full response dataset using participant-specific fitted parameters in CBM.

N = numel(templateData);
M = numel(templateData{1}.y);

modelChoicePerSubj = modelChoicePerSubj(:);
Ysim = nan(N, M);

for s = 1:N
    m = modelChoicePerSubj(s);
    par = get_subject_parameters(cbm, m, s);

    for i = 1:M
        p = compute_predicted_bin_probs(m, par, i);
        Ysim(s,i) = draw_categorical(p);
    end
end

sim_data = cell(N, 1);
for s = 1:N
    sim_data{s} = struct('y', Ysim(s,:));
end
end

function par = get_subject_parameters(cbm, modelIdx, subjIdx)
P = cbm.output.parameters{modelIdx};

if isempty(P)
    par = [];
else
    par = P(subjIdx,:);
end

end

function p = compute_predicted_bin_probs(modelIdx, params, trialIdx)
%COMPUTE_PREDICTED_BIN_PROBS  Discrete probabilities over six Likert bins.
%
%   P = COMPUTE_PREDICTED_BIN_PROBS(MODELIDX, PARAMS, TRIALIDX) returns a
%   1x6 vector of probabilities for Likert response categories 1 through 6
%   under the specified behavioural model and parameterisation.

if trialIdx < 1 || trialIdx > 16
    error('compute_predicted_bin_probs: trialIdx must be in 1..16');
end

switch modelIdx
    case 1  % mentalizing
        % Params: [log_ce]
        ce = exp(params(1));
        % Weight list for experiment 3 conditions
        w_list = [0.0000000, 0.0000000, 1.0000000, 1.0000000, ...
                  0.2735507, 0.7765487, 0.4667553, 0.5992063, ...
                  0.0000000, 0.0000000, 1.0000000, 1.0000000, ...
                  0.2735507, 0.7765487, 0.4667553, 0.5992063];
        w = w_list(trialIdx);
        alpha = 1 + w  * ce;
        beta  = 1 + (1 - w) * ce;
        % Cumulative distribution at bin edges 0,1/6,...,1
        F = betainc((0:6)/6, alpha, beta);
        p = diff(F);

    case 2  % level-0 vt_model
        % Params: [log_ce, log_n]
        ce = exp(params(1));
        n  = exp(params(2));
        denom = 2 + 2*n;
        w_list2 = [1/denom, (1+2*n)/denom, 0.5, 0.5, ...
                   1/denom, (1+2*n)/denom, 1/denom, (1+2*n)/denom, ...
                   1/denom, (1+2*n)/denom, 0.5, 0.5, ...
                   1/denom, (1+2*n)/denom, 1/denom, (1+2*n)/denom];
        w_list1 = [(1+2*n)/denom, 1/denom, 0.5, 0.5, ...
                   (1+2*n)/denom, 1/denom, (1+2*n)/denom, 1/denom, ...
                   (1+2*n)/denom, 1/denom, 0.5, 0.5, ...
                   (1+2*n)/denom, 1/denom, (1+2*n)/denom, 1/denom];
        w1 = w_list1(trialIdx);
        w2 = w_list2(trialIdx);
        alpha = 1 + w2 * ce;
        beta  = 1 + w1 * ce;
        F = betainc((0:6)/6, alpha, beta);
        p = diff(F);

    case 3  % level-1 dw_model
        % Params: [log_ce, log_n, log_e]
        ce = exp(params(1));
        n  = exp(params(2));
        e  = exp(params(3));
        denom_n  = 2 + 2*n;
        denom_ne = 2 + 2*n*e;
        w_list2 = [1/denom_ne, (1+2*n*e)/denom_ne, 0.5, 0.5, ...
                   1/denom_n,  (1+2*n)/denom_n,    1/denom_n, (1+2*n)/denom_n, ...
                   1/denom_ne, (1+2*n*e)/denom_ne, 0.5, 0.5, ...
                   1/denom_n,  (1+2*n)/denom_n,    1/denom_n, (1+2*n)/denom_n];
        w_list1 = [(1+2*n*e)/denom_ne, 1/denom_ne, 0.5, 0.5, ...
                   (1+2*n)/denom_n,    1/denom_n,  (1+2*n)/denom_n, 1/denom_n, ...
                   (1+2*n*e)/denom_ne, 1/denom_ne, 0.5, 0.5, ...
                   (1+2*n)/denom_n,    1/denom_n,  (1+2*n)/denom_n, 1/denom_n];
        w1 = w_list1(trialIdx);
        w2 = w_list2(trialIdx);
        alpha = 1 + w2 * ce;
        beta  = 1 + w1 * ce;
        F = betainc((0:6)/6, alpha, beta);
        p = diff(F);

    case 4  % level-2 kw_model
        % Params: [log_ce, log_n]
        ce = exp(params(1));
        n  = exp(params(2));
        % Colour constants for conditions
        eargain  = 1.675676;
        footgain = 1.027027;
        % Construct weight list (16 entries)
        w_list = [ ...
            1 / (2 + (eargain + footgain) * n), ...
            (1 + (eargain + footgain) * n) / (2 + (eargain + footgain) * n), ...
            (1 + footgain * n) / (2 + (eargain + footgain) * n), ...
            (1 + eargain * n) / (2 + (eargain + footgain) * n), ...
            1 / (2 + eargain * 2 * n), ...
            (1 + eargain * 2 * n) / (2 + eargain * 2 * n), ...
            1 / (2 + footgain * 2 * n), ...
            (1 + footgain * 2 * n) / (2 + footgain * 2 * n), ...
            1 / (2 + (eargain + footgain) * n), ...
            (1 + (eargain + footgain) * n) / (2 + (eargain + footgain) * n), ...
            (1 + eargain * n) / (2 + (eargain + footgain) * n), ...
            (1 + footgain * n) / (2 + (eargain + footgain) * n), ...
            1 / (2 + eargain * 2 * n), ...
            (1 + eargain * 2 * n) / (2 + eargain * 2 * n), ...
            1 / (2 + footgain * 2 * n), ...
            (1 + footgain * 2 * n) / (2 + footgain * 2 * n) ...
            ];
        w1 = w_list(trialIdx);
        w2 = 1 - w1;
        alpha = 1 + w1 * ce;
        beta  = 1 + w2 * ce;
        F = betainc((0:6)/6, alpha, beta);
        p = diff(F);

    case 5  % baseline uniform
        p = ones(1,6) / 6;

    otherwise
        error('compute_predicted_bin_probs: unknown modelIdx %d', modelIdx);
end

% Protect against numerical issues
p = max(p, realmin);
% Normalise to ensure sum(p) = 1
p = p / sum(p);
end

function k = draw_categorical(p)
c = cumsum(p(:));
u = rand();
k = find(u <= c, 1, 'first');
if isempty(k)
    k = numel(p);
end
end

function r2 = pooled_r2(Yobs, Ysim, valid_mask)
% Calculate r2 by pooling all trials across all participants.
valid = valid_mask & isfinite(Ysim);
obs = Yobs(valid);
sim = Ysim(valid);
obs = obs(:);
sim = sim(:);

if numel(obs) < 2 || std(obs) == 0 || std(sim) == 0
    r2 = NaN;
else
    C = corrcoef(obs, sim);
    r2 = C(1, 2)^2;
end
end

function r2 = subject_r2(Yobs, Ysim, valid_mask)
% calculate r2 separately for each participant, then average across participants.
n_participants = size(Yobs, 1);
r2 = nan(n_participants, 1);
for s = 1:n_participants
    valid = valid_mask(s, :) & isfinite(Ysim(s, :));
    obs_s = Yobs(s, valid);
    sim_s = Ysim(s, valid);
    obs_s = obs_s(:);
    sim_s = sim_s(:);

    if numel(obs_s) >= 2 && std(obs_s) > 0 && std(sim_s) > 0
        C = corrcoef(obs_s, sim_s);
        r2(s) = C(1, 2)^2;
    end
end
end
