function ll = level1_model_exp3(parameter, subj)
%LEVEL1_MODEL_EXP3  Log‑likelihood for the diversity heuristic (level‑1) model (Experiment 3).
%
%   LL = LEVEL1_MODEL_EXP3(PARAMETER, SUBJ) returns the log‑likelihood of
%   the observed responses under the diversity heuristic model for
%   Experiment 3.  This model has three free parameters: the log
%   concentration (log_ce), the log noise term (log_n) and the log
%   elasticity (log_e).  The model uses two sets of denominators
%   depending on n and n*e to build weight lists for the 16 trials.
%
%   Inputs:
%     PARAMETER – [log_ce, log_n, log_e].
%     SUBJ      – struct with field ``y`` containing a 1×16 vector of
%                 Likert responses (1…6) or NaN.
%
%   Output:
%     LL        – scalar log‑likelihood of the data under the model.
%
%   See also MENTALIZING_MODEL_EXP3, LEVEL0_MODEL_EXP3, LEVEL2_MODEL_EXP3,
%            BASELINE_MODEL_EXP3.

% Exponentiate parameters
ce = exp(parameter(1));
n  = exp(parameter(2));
e  = exp(parameter(3));

responses = subj.y;

% Denominators for n and n*e
denom_n  = 2 + 2 * n;
denom_ne = 2 + 2 * n * e;

% Define weight lists for 16 conditions.  The pattern alternates between
% conditions that depend on n*e and n, reflecting the diversity heuristic.
w_list2 = [1/denom_ne, (1+2*n*e)/denom_ne, 0.5, 0.5, ...
           1/denom_n,  (1+2*n)/denom_n,    1/denom_n, (1+2*n)/denom_n, ...
           1/denom_ne, (1+2*n*e)/denom_ne, 0.5, 0.5, ...
           1/denom_n,  (1+2*n)/denom_n,    1/denom_n, (1+2*n)/denom_n];
w_list1 = [(1+2*n*e)/denom_ne, 1/denom_ne, 0.5, 0.5, ...
           (1+2*n)/denom_n,    1/denom_n,  (1+2*n)/denom_n, 1/denom_n, ...
           (1+2*n*e)/denom_ne, 1/denom_ne, 0.5, 0.5, ...
           (1+2*n)/denom_n,    1/denom_n,  (1+2*n)/denom_n, 1/denom_n];

% Initialise log‑likelihood
ll = 0;

for i = 1:numel(responses)
    y = responses(i);
    if isnan(y) || y < 1 || y > 6
        continue;
    end
    w2 = w_list2(i);
    w1 = w_list1(i);
    alpha = 1 + w2 * ce;
    beta  = 1 + w1 * ce;
    p = diff(betainc((0:6)/6, alpha, beta));
    p = max(p, realmin);
    p = p / sum(p);
    ll = ll + log(p(round(y)));
end
end