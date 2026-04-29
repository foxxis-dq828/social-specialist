function [mean_r_squared, r_squared] = get_mean_r_squared_bmi_exp3(nboot, data, cbm, modelChoicePerSubj)
%GET_MEAN_R_SQUARED_BMI_EXP3  Pooled R² between observed and simulated data (Experiment 3).
%
%   [MEAN_R_SQUARED, R_SQUARED] = GET_MEAN_R_SQUARED_BMI_EXP3(NBOOT, DATA, CBM, MODELCHOICEPERSUBJ)
%   repeatedly draws synthetic datasets from the posterior predictive
%   distribution of each subject's selected model (for Experiment 3) and
%   computes the squared Pearson correlation coefficient between the
%   observed responses and the simulated responses, pooling over all
%   subjects and trials.  The function returns the mean pooled R² across
%   bootstrap iterations (MEAN_R_SQUARED) and the vector of bootstrapped
%   R² values (R_SQUARED).
%
%   Inputs:
%     NBOOT             – number of bootstrap simulations to draw.
%     DATA              – cell array of subject structs with field ``y``.
%     CBM               – output structure from cbm_hbi containing group
%                          and individual parameter estimates.
%     MODELCHOICEPERSUBJ– N×1 vector indicating which model index applies
%                          to each subject.
%
%   See also GET_INDIVIDUAL_MEAN_R_SQUARED_BMI_EXP3.

disp('Calculating mean r‑squared iteratively (pooled across subjects) for Experiment 3');

N = numel(data);
M = numel(data{1}.y);

% Build observed matrix N×M
Yobs = nan(N, M);
for s = 1:N
    y = data{s}.y(:)';
    if numel(y) ~= M
        error('Subject %d has %d trials; expected %d. Ensure equal‑length y.', s, numel(y), M);
    end
    Yobs(s,:) = y;
end

% Mask valid observations
validMask = ~isnan(Yobs) & (Yobs >= 1) & (Yobs <= 6);

r_squared = nan(nboot,1);

for boot = 1:nboot
    % Print progress every 10 iterations
    if mod(boot, 10) == 0
        disp(['Bootstrap iteration ' num2str(boot)]);
    end

    Ysim = simulate_dataset_from_hbi_local(cbm, modelChoicePerSubj, N, M);

    obs = Yobs(validMask);
    sim = Ysim(validMask);
    obs = obs(:);
    sim = sim(:);

    if numel(obs) < 2 || std(obs) == 0 || std(sim) == 0
        r_squared(boot) = NaN;
    else
        r = corr(obs, sim);
        r_squared(boot) = r^2;
    end
end

mean_r_squared = mean(r_squared, 'omitnan');
end

%% ============================ Local helpers ============================

function Ysim = simulate_dataset_from_hbi_local(cbm, modelChoicePerSubj, N, M)
%SIMULATE_DATASET_FROM_HBI_LOCAL  Generate synthetic responses from HBI posterior (Experiment 3).
Ysim = nan(N, M);

for s = 1:N
    m = modelChoicePerSubj(s);

    % Extract individual parameters.  cbm.output.parameters{m} may be a cell or matrix.
    P = cbm.output.parameters{m};
    if isempty(P)
        par = [];
    elseif iscell(P)
        par = P{s};
        if isrow(par)
            par = par(:)';
        end
    else
        par = P(s, :);
    end

    for i = 1:M
        p = model_bin_probs_local(m, par, i);
        Ysim(s, i) = draw_categorical_local(p);
    end
end
end

function p = model_bin_probs_local(modelIdx, par, itemIdx)
%MODEL_BIN_PROBS_LOCAL  Local predicted probabilities for Experiment 3.

switch modelIdx
    case 1  % Mentalizing
        ce = exp(clamp_local(par(1)));
        w_list = [0.0000000, 0.0000000, 1.0000000, 1.0000000, ...
                  0.2735507, 0.7765487, 0.4667553, 0.5992063, ...
                  0.0000000, 0.0000000, 1.0000000, 1.0000000, ...
                  0.2735507, 0.7765487, 0.4667553, 0.5992063];
        w = w_list(itemIdx);
        alpha = 1 + w * ce;
        beta  = 1 + (1 - w) * ce;
        cdfVals = betacdf((0:6)/6, alpha, beta);
        p = diff(cdfVals);
    case 2  % Level‑0
        ce = exp(clamp_local(par(1)));
        n  = exp(clamp_local(par(2)));
        denom = 2 + 2 * n;
        w_list2 = [1/denom, (1+2*n)/denom, 0.5, 0.5, ...
                   1/denom, (1+2*n)/denom, 1/denom, (1+2*n)/denom, ...
                   1/denom, (1+2*n)/denom, 0.5, 0.5, ...
                   1/denom, (1+2*n)/denom, 1/denom, (1+2*n)/denom];
        w_list1 = [(1+2*n)/denom, 1/denom, 0.5, 0.5, ...
                   (1+2*n)/denom, 1/denom, (1+2*n)/denom, 1/denom, ...
                   (1+2*n)/denom, 1/denom, 0.5, 0.5, ...
                   (1+2*n)/denom, 1/denom, (1+2*n)/denom, 1/denom];
        w2 = w_list2(itemIdx);
        w1 = w_list1(itemIdx);
        alpha = 1 + w2 * ce;
        beta  = 1 + w1 * ce;
        F = betainc((0:6)/6, alpha, beta);
        p = diff(F);
    case 3  % Diversity heuristic
        ce = exp(clamp_local(par(1)));
        n  = exp(clamp_local(par(2)));
        e  = exp(clamp_local(par(3)));
        denom_n  = 2 + 2 * n;
        denom_ne = 2 + 2 * n * e;
        w_list2 = [1/denom_ne, (1+2*n*e)/denom_ne, 0.5, 0.5, ...
                   1/denom_n,  (1+2*n)/denom_n,    1/denom_n, (1+2*n)/denom_n, ...
                   1/denom_ne, (1+2*n*e)/denom_ne, 0.5, 0.5, ...
                   1/denom_n,  (1+2*n)/denom_n,    1/denom_n, (1+2*n)/denom_n];
        w_list1 = [(1+2*n*e)/denom_ne, 1/denom_ne, 0.5, 0.5, ...
                   (1+2*n)/denom_n,    1/denom_n,  (1+2*n)/denom_n, 1/denom_n, ...
                   (1+2*n*e)/denom_ne, 1/denom_ne, 0.5, 0.5, ...
                   (1+2*n)/denom_n,    1/denom_n,  (1+2*n)/denom_n, 1/denom_n];
        w2 = w_list2(itemIdx);
        w1 = w_list1(itemIdx);
        alpha = 1 + w2 * ce;
        beta  = 1 + w1 * ce;
        F = betainc((0:6)/6, alpha, beta);
        p = diff(F);
    case 4  % Colour model
        ce = exp(clamp_local(par(1)));
        n  = exp(clamp_local(par(2)));
        eargain  = 1.675676;
        footgain = 1.027027;
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
        w1 = w_list(itemIdx);
        w2 = 1 - w1;
        alpha = 1 + w1 * ce;
        beta  = 1 + w2 * ce;
        F = betainc((0:6)/6, alpha, beta);
        p = diff(F);
    case 5  % Baseline symmetric model
        ce = exp(clamp_local(par(1)));
        w  = 0.5;
        alpha = 1 + w * ce;
        beta  = 1 + (1 - w) * ce;
        F = betainc((0:6)/6, alpha, beta);
        p = diff(F);
    otherwise
        error('Unknown modelIdx: %d', modelIdx);
end

% Protect and normalise
p = max(p, realmin);
p = p / sum(p);
end

function k = draw_categorical_local(p)
%DRAW_CATEGORICAL_LOCAL  Draw a sample from a categorical distribution.
c = cumsum(p(:));
u = rand();
k = find(u <= c, 1, 'first');
if isempty(k)
    k = numel(p);
end
end

function z = clamp_local(z)
%CLAMP_LOCAL  Clamp a value into [-20, 20] to avoid numerical issues.
z = min(max(z, -20), 20);
end