function [mean_r_squared, r_squared] = get_mean_r_squared_bmi_exp2(nboot, data, cbm, modelChoicePerSubj)
%GET_MEAN_R_SQUARED_BMI_EXP2  Bootstrapped pooled R² across all subjects (Experiment 2).
%
%   [MEAN_R_SQUARED, R_SQUARED] = GET_MEAN_R_SQUARED_BMI_EXP2(NBOOT, DATA, CBM, MODELCHOICEPERSUBJ)
%   repeatedly simulates synthetic response datasets from the posterior
%   predictive distribution of each subject's selected model (for Experiment 2)
%   and computes the squared Pearson correlation coefficient (R²) between
%   the observed data and the simulated data, pooling all subjects together.
%   The function returns the mean pooled R² across bootstrap iterations
%   (MEAN_R_SQUARED) and the vector of bootstrap R² values (R_SQUARED).
%
%   Inputs:
%     NBOOT             – number of bootstrap simulations to draw.
%     DATA              – cell array of subject structs with field ``y``.
%     CBM               – output structure from cbm_hbi containing group and
%                          individual parameter estimates.
%     MODELCHOICEPERSUBJ– N×1 vector indicating which model index applies to each
%                          subject.
%
%   See also GET_INDIVIDUAL_MEAN_R_SQUARED_BMI_EXP2.

disp('Calculating mean r‑squared iteratively (pooled across subjects)')

N = numel(data);
M = numel(data{1}.y);

% Build observed matrix N×M
Yobs = nan(N,M);
for s = 1:N
    y = data{s}.y(:)';
    if numel(y) ~= M
        error('Subject %d has %d trials; expected %d. Ensure equal-length y.', s, numel(y), M);
    end
    Yobs(s,:) = y;
end

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

    if numel(obs) < 2 || std(obs)==0 || std(sim)==0
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
Ysim = nan(N,M);

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
        par = P(s,:);
    end

    for i = 1:M
        p = model_bin_probs_local(m, par, i);
        Ysim(s,i) = draw_categorical_local(p);
    end
end
end

function p = model_bin_probs_local(modelIdx, par, itemIdx)

% Local predicted probability function specific to Experiment 2
switch modelIdx
    case 1
        ce = exp(clamp_local(par(1)));
        i = itemIdx;
        w_list = [0.6897720, 0.8861646, 0.3220230, 0.6803653, 0.5384914, 0.5425532, 0.2733100, 0.3704918];
        w     = w_list(i);
        alpha = 1 + w * ce;
        beta  = 1 + (1 - w) * ce;
        cdfVals = betainc((0:6)/6, alpha, beta);
        p = diff(cdfVals);

    case 2
        ce = exp(clamp_local(par(1)));
        n  = exp(clamp_local(par(2)));
        i = itemIdx;
        denom = 2 + 2*n;
        w_list2 = [(1+2*n)/denom, (1+2*n)/denom, 1/denom, 0.5, (1+2*n)/denom, (1+2*n)/denom, 1/denom, 0.5];
        w_list1 = [1/denom, 1/denom, (1+2*n)/denom, 0.5, 1/denom, 1/denom, (1+2*n)/denom, 0.5];
        w1 = w_list1(i);
        w2 = w_list2(i);
        alpha = 1 + w2 * ce;
        beta  = 1 + w1 * ce;
        F = betainc((0:6)/6, alpha, beta);
        p = diff(F);

    case 3
        ce = exp(clamp_local(par(1)));
        n  = exp(clamp_local(par(2)));
        e  = exp(clamp_local(par(3)));
        i = itemIdx;
        denom_n  = 2 + 2*n;
        denom_ne = 2 + 2*n*e;
        w_list2 = [(1+2*n)/denom_n, (1+2*n*e)/denom_ne, 1/denom_n, 0.5, (1+2*n)/denom_n, (1+2*n*e)/denom_ne, 1/denom_n, 0.5];
        w_list1 = [1/denom_n, 1/denom_ne, (1+2*n)/denom_n, 0.5, 1/denom_n, 1/denom_ne, (1+2*n)/denom_n, 0.5];
        w1 = w_list1(i);
        w2 = w_list2(i);
        alpha = 1 + w2 * ce;
        beta  = 1 + w1 * ce;
        F = betainc((0:6)/6, alpha, beta);
        p = diff(F);

    case 4
        ce = exp(clamp_local(par(1)));
        n  = exp(clamp_local(par(2)));
        i = itemIdx;
        d_blue   = [1.604278, 1.604278, 1.410604, 1.410604, 1.040060, 1.040060, 1.463840, 1.463840];
        d_yellow = [1.604278, 1.604278, 1.410604, 1.410604, 1.040060, 1.040060, 1.463840, 1.463840];
        w_list = zeros(1, 8);
        w_list(1) = (1 + d_blue(1)*2*n) / (2 + d_blue(1)*2*n);
        w_list(2) = (1 + (d_blue(2) + d_yellow(2))*n) / (2 + (d_blue(2) + d_yellow(2))*n);
        w_list(3) = 1 / (2 + d_yellow(3)*2*n);
        w_list(4) = (1 + d_blue(4)*n) / (2 + (d_yellow(4) + d_blue(4))*n);
        w_list(5) = (1 + d_blue(5)*2*n) / (2 + d_blue(5)*2*n);
        w_list(6) = (1 + (d_blue(6) + d_yellow(6))*n) / (2 + (d_blue(6) + d_yellow(6))*n);
        w_list(7) = 1 / (2 + d_blue(7)*2*n);
        w_list(8) = (1 + d_yellow(8)*n) / (2 + (d_yellow(8) + d_blue(8))*n);
        w1 = w_list(i);
        w2 = 1 - w1;
        alpha = 1 + w1 * ce;
        beta  = 1 + w2 * ce;
        F = betainc((0:6)/6, alpha, beta);
        p = diff(F);

    case 5
        p  = ones(1,6) / 6;

    otherwise
        error('Unknown modelIdx: %d', modelIdx);
end

% Protect and normalise
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
