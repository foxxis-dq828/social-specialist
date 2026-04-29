function ll = level1_model_exp2(parameter, subj)
%LEVEL1_MODEL_EXP2  Log‑likelihood for the diversity heuristic (level‑1) model (Experiment 2).
%
%   LL = LEVEL1_MODEL_EXP2(PARAMETER, SUBJ) returns the log‑likelihood of
%   a subject's responses under the diversity heuristic model.  PARAMETER
%   is [log_ce, log_n, log_e], representing the log of the concentration
%   parameter ce, the noisiness parameter n and the elasticity parameter e.
%
%   The model uses two weight lists that depend on n and n*e to compute
%   beta distribution parameters for each trial.
%
%   Inputs:
%     PARAMETER – [log_ce, log_n, log_e].
%     SUBJ      – struct with field y.
%
%   Output:
%     LL – scalar log‑likelihood.
%
%   See also MENTALIZING_MODEL_EXP2, LEVEL0_MODEL_EXP2, LEVEL2_MODEL_EXP2.

ce = exp(parameter(1));
n  = exp(parameter(2));
e  = exp(parameter(3));
responses = subj.y;

ll = 0;
denom_n  = 2 + 2*n;
denom_ne = 2 + 2*n*e;

w_list2 = [(1+2*n)/denom_n, (1+2*n*e)/denom_ne, 1/denom_n, 0.5, ...
           (1+2*n)/denom_n, (1+2*n*e)/denom_ne, 1/denom_n, 0.5];
w_list1 = [1/denom_n, 1/denom_ne, (1+2*n)/denom_n, 0.5, ...
           1/denom_n, 1/denom_ne, (1+2*n)/denom_n, 0.5];

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
