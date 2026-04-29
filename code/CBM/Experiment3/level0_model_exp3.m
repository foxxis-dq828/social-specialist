function ll = level0_model_exp3(parameter, subj)
%LEVEL0_MODEL_EXP3  Log‑likelihood for the level‑0 model (Experiment 3).
%
%   LL = LEVEL0_MODEL_EXP3(PARAMETER, SUBJ) computes the log‑likelihood of
%   a subject's responses under the level‑0 model used in Experiment 3.
%   This model has two free parameters: the log concentration (log_ce) and
%   the log noise term (log_n).  The weights determining the means of the
%   Beta distributions repeat a pattern of eight conditions twice to
%   accommodate the 16 conditions in Experiment 3.  Weight lists w_list1
%   and w_list2 correspond to the Beta parameters beta and alpha,
%   respectively.
%
%   Inputs:
%     PARAMETER – [log_ce, log_n].
%     SUBJ      – struct with field ``y`` containing a 1×16 vector of
%                 Likert responses or NaN.
%
%   Output:
%     LL        – scalar log‑likelihood of the data under the model.
%
%   See also MENTALIZING_MODEL_EXP3, LEVEL1_MODEL_EXP3, LEVEL2_MODEL_EXP3,
%            BASELINE_MODEL_EXP3.

% Unpack and exponentiate parameters
ce = exp(parameter(1));
n  = exp(parameter(2));

responses = subj.y;

% Precompute denominators for the weight lists
denom = 2 + 2 * n;

% Define weight lists for 16 conditions.  The pattern for eight conditions
% is repeated to cover all 16 trials.  w_list2 affects the alpha (first
% Beta parameter) and w_list1 affects the beta (second parameter).
w_list2 = [1/denom, (1+2*n)/denom, 0.5, 0.5, ...
           1/denom, (1+2*n)/denom, 1/denom, (1+2*n)/denom, ...
           1/denom, (1+2*n)/denom, 0.5, 0.5, ...
           1/denom, (1+2*n)/denom, 1/denom, (1+2*n)/denom];
w_list1 = [(1+2*n)/denom, 1/denom, 0.5, 0.5, ...
           (1+2*n)/denom, 1/denom, (1+2*n)/denom, 1/denom, ...
           (1+2*n)/denom, 1/denom, 0.5, 0.5, ...
           (1+2*n)/denom, 1/denom, (1+2*n)/denom, 1/denom];

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