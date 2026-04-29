function p = compute_predicted_bin_probs_exp2(modelIdx, params, trialIdx)
%COMPUTE_PREDICTED_BIN_PROBS_EXP2  Discrete probabilities over six Likert bins (Experiment 2).
%
%   P = COMPUTE_PREDICTED_BIN_PROBS_EXP2(MODELIDX, PARAMS, TRIALIDX) returns
%   a 1×6 vector of probabilities for Likert response categories 1 through 6
%   under the specified behavioural model and parameterisation for Experiment 2.
%   PARAMS is supplied in log space and TRIALIDX is the experimental condition index (1…8).
%
%   MODELIDX:
%     1 – mentalizing model (one parameter).
%     2 – level‑0 model (two parameters).
%     3 – diversity heuristic / level‑1 model (three parameters).
%     4 – colour / level‑2 model (two parameters).
%     5 – baseline uniform model.
%
%   The function normalises probabilities and protects against underflow.
%
%   See also COMPUTE_PREDICTED_MEAN_EXP2.

if trialIdx < 1 || trialIdx > 8
    error('compute_predicted_bin_probs_exp2: trialIdx must be in 1..8');
end

switch modelIdx
    case 1
        ce = exp(params(1));
        w_list = [0.6897720, 0.8861646, 0.3220230, 0.6803653, 0.5384914, 0.5425532, 0.2733100, 0.3704918];
        w = w_list(trialIdx);
        alpha = 1 + w * ce;
        beta  = 1 + (1 - w) * ce;
        F = betainc((0:6)/6, alpha, beta);
        p = diff(F);
    case 2
        ce = exp(params(1));
        n  = exp(params(2));
        denom = 2 + 2*n;
        w_list2 = [(1+2*n)/denom, (1+2*n)/denom, 1/denom, 0.5, (1+2*n)/denom, (1+2*n)/denom, 1/denom, 0.5];
        w_list1 = [1/denom, 1/denom, (1+2*n)/denom, 0.5, 1/denom, 1/denom, (1+2*n)/denom, 0.5];
        w2 = w_list2(trialIdx);
        w1 = w_list1(trialIdx);
        alpha = 1 + w2 * ce;
        beta  = 1 + w1 * ce;
        F = betainc((0:6)/6, alpha, beta);
        p = diff(F);
    case 3
        ce = exp(params(1));
        n  = exp(params(2));
        e  = exp(params(3));
        denom_n  = 2 + 2*n;
        denom_ne = 2 + 2*n*e;
        w_list2 = [(1+2*n)/denom_n, (1+2*n*e)/denom_ne, 1/denom_n, 0.5, (1+2*n)/denom_n, (1+2*n*e)/denom_ne, 1/denom_n, 0.5];
        w_list1 = [1/denom_n, 1/denom_ne, (1+2*n)/denom_n, 0.5, 1/denom_n, 1/denom_ne, (1+2*n)/denom_n, 0.5];
        w2 = w_list2(trialIdx);
        w1 = w_list1(trialIdx);
        alpha = 1 + w2 * ce;
        beta  = 1 + w1 * ce;
        F = betainc((0:6)/6, alpha, beta);
        p = diff(F);
    case 4
        ce = exp(params(1));
        n  = exp(params(2));
        d_blue   = [1.604278, 1.604278, 1.410604, 1.410604, 1.040060, 1.040060, 1.463840, 1.463840];
        d_yellow = [1.604278, 1.604278, 1.410604, 1.410604, 1.040060, 1.040060, 1.463840, 1.463840];
        w_list = zeros(1,8);
        w_list(1) = (1 + d_blue(1)*2*n) / (2 + d_blue(1)*2*n);
        w_list(2) = (1 + (d_blue(2) + d_yellow(2))*n) / (2 + (d_blue(2) + d_yellow(2))*n);
        w_list(3) = 1 / (2 + d_yellow(3)*2*n);
        w_list(4) = (1 + d_blue(4)*n) / (2 + (d_yellow(4) + d_blue(4))*n);
        w_list(5) = (1 + d_blue(5)*2*n) / (2 + d_blue(5)*2*n);
        w_list(6) = (1 + (d_blue(6) + d_yellow(6))*n) / (2 + (d_blue(6) + d_yellow(6))*n);
        w_list(7) = 1 / (2 + d_blue(7)*2*n);
        w_list(8) = (1 + d_yellow(8)*n) / (2 + (d_yellow(8) + d_blue(8))*n);
        w1 = w_list(trialIdx);
        w2 = 1 - w1;
        alpha = 1 + w1 * ce;
        beta  = 1 + w2 * ce;
        F = betainc((0:6)/6, alpha, beta);
        p = diff(F);
    case 5
        p = ones(1,6) / 6;
    otherwise
        error('compute_predicted_bin_probs_exp2: unknown modelIdx %d', modelIdx);
end

p = max(p, realmin);
p = p / sum(p);
end
