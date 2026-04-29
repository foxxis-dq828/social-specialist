function [mean_r2_indiv, r2_indiv_boot] = get_individual_mean_r_squared_bmi_exp3(nboot, data, cbm, modelChoicePerSubj)
%GET_INDIVIDUAL_MEAN_R_SQUARED_BMI_EXP3  Subject‑wise R² values for Experiment 3.
%
%   [MEAN_R2_INDIV, R2_INDIV_BOOT] = GET_INDIVIDUAL_MEAN_R_SQUARED_BMI_EXP3(NBOOT, DATA, CBM, MODELCHOICEPERSUBJ)
%   draws NBOOT synthetic datasets from the posterior predictive
%   distribution of each subject's chosen model (for Experiment 3) and
%   computes the squared Pearson correlation coefficient between observed
%   and simulated responses on a per‑subject basis.  The function returns
%   the average of these individual R² values across subjects and the
%   vector of bootstrapped means across iterations.
%
%   Inputs:
%     NBOOT             – number of bootstrap simulations.
%     DATA              – cell array of subject structs with field ``y``.
%     CBM               – CBM output structure from cbm_hbi.
%     MODELCHOICEPERSUBJ– N×1 vector indicating which model applies to each
%                          subject.
%
%   Outputs:
%     MEAN_R2_INDIV     – mean of individual R² values across subjects and
%                          bootstrap samples.
%     R2_INDIV_BOOT     – vector of mean R² values for each bootstrap iteration.
%
%   See also GET_MEAN_R_SQUARED_BMI_EXP3.

disp('Calculating mean r‑squared iteratively (individual level) for Experiment 3');

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

validMask = ~isnan(Yobs) & (Yobs >= 1) & (Yobs <= 6);

r2_indiv_boot = nan(nboot, 1);
r2_subj_boot  = nan(N, nboot);

for boot = 1:nboot
    if mod(boot, 10) == 0
        disp(['Bootstrap iteration ' num2str(boot)]);
    end
    Ysim = simulate_dataset_from_hbi_local(cbm, modelChoicePerSubj, N, M);
    % Compute per‑subject R²
    for s = 1:N
        v = validMask(s, :);
        obs_s = Yobs(s, v);
        sim_s = Ysim(s, v);
        if numel(obs_s) < 2 || std(obs_s) == 0 || std(sim_s) == 0
            r2_subj_boot(s, boot) = NaN;
        else
            r = corr(obs_s(:), sim_s(:));
            r2_subj_boot(s, boot) = r^2;
        end
    end
    r2_indiv_boot(boot) = mean(r2_subj_boot(:, boot), 'omitnan');
end

mean_r2_indiv = mean(r2_indiv_boot, 'omitnan');
end

%% ============================ Local helpers ============================

function Ysim = simulate_dataset_from_hbi_local(cbm, modelChoicePerSubj, N, M)
Ysim = nan(N, M);
for s = 1:N
    m = modelChoicePerSubj(s);
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
switch modelIdx
    case 1
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
    case 2
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
    case 3
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
    case 4
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
    case 5
        ce = exp(clamp_local(par(1)));
        w  = 0.5;
        alpha = 1 + w * ce;
        beta  = 1 + (1 - w) * ce;
        F = betainc((0:6)/6, alpha, beta);
        p = diff(F);
    otherwise
        error('Unknown modelIdx: %d', modelIdx);
end
p = max(p, realmin);
p = p / sum(p);
end

function k = draw_categorical_local(p)
c = cumsum(p(:));
u = rand();
k = find(u <= c, 1, 'first');
if isempty(k)
    k = numel(p);
end
end

function z = clamp_local(z)
z = min(max(z, -20), 20);
end