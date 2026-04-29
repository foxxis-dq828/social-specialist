function p = compute_predicted_bin_probs_exp3(modelIdx, params, trialIdx)
%COMPUTE_PREDICTED_BIN_PROBS_EXP3  Predicted response probabilities for Experiment 3.
%
%   P = COMPUTE_PREDICTED_BIN_PROBS_EXP3(MODELIDX, PARAMS, TRIALIDX) returns
%   a 1×6 vector of probabilities for Likert response categories 1–6
%   under the specified behavioural model and parameterisation for
%   Experiment 3.  PARAMS is supplied in log space.  TRIALIDX must be an
%   integer in 1…16.
%
%   MODELIDX:
%     1 – mentalizing model (one parameter).
%     2 – level‑0 model (two parameters).
%     3 – diversity heuristic / level‑1 model (three parameters).
%     4 – colour / level‑2 model (two parameters).
%     5 – baseline symmetric model (one parameter).
%
%   The function normalises probabilities and protects against underflow.
%
%   See also COMPUTE_PREDICTED_MEAN_EXP3.

% Validate trial index
if trialIdx < 1 || trialIdx > 16
    error('compute_predicted_bin_probs_exp3: trialIdx must be in 1..16');
end

switch modelIdx
    case 1  % Mentalizing model
        ce = exp(params(1));
        w_list = [0.0000000, 0.0000000, 1.0000000, 1.0000000, ...
                  0.2735507, 0.7765487, 0.4667553, 0.5992063, ...
                  0.0000000, 0.0000000, 1.0000000, 1.0000000, ...
                  0.2735507, 0.7765487, 0.4667553, 0.5992063];
        w = w_list(trialIdx);
        alpha = 1 + w * ce;
        beta  = 1 + (1 - w) * ce;
        F = betainc((0:6)/6, alpha, beta);
        p = diff(F);
    case 2  % Level‑0 model
        ce = exp(params(1));
        n  = exp(params(2));
        denom = 2 + 2 * n;
        w_list2 = [1/denom, (1+2*n)/denom, 0.5, 0.5, ...
                   1/denom, (1+2*n)/denom, 1/denom, (1+2*n)/denom, ...
                   1/denom, (1+2*n)/denom, 0.5, 0.5, ...
                   1/denom, (1+2*n)/denom, 1/denom, (1+2*n)/denom];
        w_list1 = [(1+2*n)/denom, 1/denom, 0.5, 0.5, ...
                   (1+2*n)/denom, 1/denom, (1+2*n)/denom, 1/denom, ...
                   (1+2*n)/denom, 1/denom, 0.5, 0.5, ...
                   (1+2*n)/denom, 1/denom, (1+2*n)/denom, 1/denom];
        w2 = w_list2(trialIdx);
        w1 = w_list1(trialIdx);
        alpha = 1 + w2 * ce;
        beta  = 1 + w1 * ce;
        F = betainc((0:6)/6, alpha, beta);
        p = diff(F);
    case 3  % Diversity heuristic (level‑1) model
        ce = exp(params(1));
        n  = exp(params(2));
        e  = exp(params(3));
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
        w2 = w_list2(trialIdx);
        w1 = w_list1(trialIdx);
        alpha = 1 + w2 * ce;
        beta  = 1 + w1 * ce;
        F = betainc((0:6)/6, alpha, beta);
        p = diff(F);
    case 4  % Colour (level‑2) model
        ce = exp(params(1));
        n  = exp(params(2));
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
        w1 = w_list(trialIdx);
        w2 = 1 - w1;
        alpha = 1 + w1 * ce;
        beta  = 1 + w2 * ce;
        F = betainc((0:6)/6, alpha, beta);
        p = diff(F);
    case 5  % Baseline symmetric model
        ce = exp(params(1));
        w  = 0.5;
        alpha = 1 + w * ce;
        beta  = 1 + (1 - w) * ce;
        F = betainc((0:6)/6, alpha, beta);
        p = diff(F);
    otherwise
        error('compute_predicted_bin_probs_exp3: unknown modelIdx %d', modelIdx);
end

% Normalise and protect against underflow
p = max(p, realmin);
p = p / sum(p);
end