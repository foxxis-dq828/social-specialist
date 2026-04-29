function ll = mentalizing_model_exp2(parameter, subj)
%MENTALIZING_MODEL_EXP2  Log‑likelihood for the mentalizing model (Experiment 2).
%
%   LL = MENTALIZING_MODEL_EXP2(PARAMETER, SUBJ) returns the log‑likelihood
%   of a subject's Likert responses under the mentalizing model used in
%   Experiment 2.  PARAMETER is a 1×1 vector containing the log of the
%   concentration parameter (log_ce).  The function uses a pre‑defined
%   weight list specific to Experiment 2 to compute the beta‑binomial
%   probabilities for each of the eight experimental conditions.
%
%   Input:
%     PARAMETER – [log_ce], the log of the concentration parameter.
%     SUBJ      – struct with field y (1×8 vector of responses, possibly NaN).
%
%   Output:
%     LL        – scalar log‑likelihood.
%
%   See also LEVEL0_MODEL_EXP2, LEVEL1_MODEL_EXP2, LEVEL2_MODEL_EXP2.
%
%   Author: automatically generated to mirror the structure of Experiment 1 code.

ce = exp(parameter(1));
responses = subj.y;
% Weight list for the eight conditions in Experiment 2 (mentalizing model).
w_list = [0.6897720, 0.8861646, 0.3220230, 0.6803653, 0.5384914, 0.5425532, 0.2733100, 0.3704918];

ll = 0;
for i = 1:numel(responses)
    y = responses(i);
    if isnan(y) || y < 1 || y > 6
        continue;
    end
    w = w_list(i);
    alpha = 1 + w * ce;
    beta  = 1 + (1 - w) * ce;
    % Compute probabilities at the six Likert categories.
    p = diff(betainc((0:6)/6, alpha, beta));
    p = max(p, realmin);
    p = p / sum(p);
    ll = ll + log(p(round(y)));
end
end
