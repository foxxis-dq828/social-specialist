function ll = level0_model_exp2(parameter, subj)
%LEVEL0_MODEL_EXP2  Log‑likelihood for the level‑0 model (Experiment 2).
%
%   LL = LEVEL0_MODEL_EXP2(PARAMETER, SUBJ) returns the log‑likelihood of
%   a subject's responses under the level‑0 model.  PARAMETER is a 1×2
%   vector [log_ce, log_n] representing the logarithms of the concentration
%   parameter (ce) and the noise parameter (n).  The model defines two
%   weight lists for each of the eight experimental conditions, which
%   determine the beta distribution parameters for the predicted rating
%   distribution.
%
%   Inputs:
%     PARAMETER – [log_ce, log_n].
%     SUBJ      – struct with field y.
%
%   Output:
%     LL – scalar log‑likelihood.
%
%   See also MENTALIZING_MODEL_EXP2, LEVEL1_MODEL_EXP2, LEVEL2_MODEL_EXP2.

ce = exp(parameter(1));
n  = exp(parameter(2));
responses = subj.y;

ll = 0;
% Denominator used in weight lists (same for all conditions)
denom = 2 + 2*n;
w_list2 = [(1+2*n)/denom, (1+2*n)/denom, 1/denom, 0.5, (1+2*n)/denom, (1+2*n)/denom, 1/denom, 0.5];
w_list1 = [1/denom, 1/denom, (1+2*n)/denom, 0.5, 1/denom, 1/denom, (1+2*n)/denom, 0.5];

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
