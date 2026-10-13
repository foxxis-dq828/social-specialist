function p = compute_predicted_bin_probs(modelIdx, params, trialIdx)
%COMPUTE_PREDICTED_BIN_PROBS  Discrete probabilities over six Likert bins.
%
%   P = COMPUTE_PREDICTED_BIN_PROBS(MODELIDX, PARAMS, TRIALIDX) returns a
%   1x6 vector of probabilities for Likert response categories 1 through 6
%   under the specified behavioural model and parameterisation. The
%   parameter vector PARAMS should be supplied in the logarithmic space
%   (i.e. the same format as returned by cbm.output.parameters). The
%   function supports five model indices corresponding to the mentalizing
%   model (1), level-0 model (2), diversity heuristic/level-1 model (3),
%   colour/level-2 model (4), and the baseline uniform model (5). TRIALIDX
%   specifies which of the sixteen experimental conditions (1..16) the
%   probabilities should be computed for. Probabilities are normalised and
%   protected against underflow.

% Ensure trial index is valid
if trialIdx < 1 || trialIdx > 16
    error('compute_predicted_bin_probs: trialIdx must be in 1..16');
end

switch modelIdx
    case 1  % mentalizing
        % Params: [log_ce]
        ce = exp(params(1));
        % Weight list for experiment 3 conditions
        w_list = [0.0000000, 0.0000000, 1.0000000, 1.0000000, ...
                  0.2727273, 0.7777778, 0.4666667, 0.6000000, ...
                  0.0000000, 0.0000000, 1.0000000, 1.0000000, ...
                  0.2727273, 0.7777778, 0.4666667, 0.6000000];
        w = w_list(trialIdx);
        alpha = 1 + w  * ce;
        beta  = 1 + (1 - w) * ce;
        % Cumulative distribution at bin edges 0,1/6,...,1
        F = betainc((0:6)/6, alpha, beta);
        p = diff(F);

    case 2  % level-0
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

    case 3  % diversity heuristic / level-1
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

    case 4  % colour / level-2
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
            (1 + footgain * n) / (2 + (eargain + footgain) * n), ...
            (1 + eargain * n) / (2 + (eargain + footgain) * n), ...
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
