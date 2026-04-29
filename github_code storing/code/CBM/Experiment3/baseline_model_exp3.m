function ll = baseline_model_exp3(parameter, subj)
%BASELINE_MODEL_EXP3  Log‑likelihood for the baseline model (Experiment 3).
%
%   LL = BASELINE_MODEL_EXP3(PARAMETER, SUBJ) computes the log‑likelihood
%   of a subject's responses under a baseline model for Experiment 3.
%   This baseline assumes that the rating distribution is symmetric around
%   the middle of the scale and is governed by a single concentration
%   parameter ce = exp(log_ce).  Specifically, a Beta distribution with
%   equal weights (w = 0.5) is used to generate probabilities for the six
%   Likert categories.  Larger ce values make the distribution more
%   concentrated around the centre, whereas ce → 0 yields a uniform
%   distribution.
%
%   Inputs:
%     PARAMETER – [log_ce], the logarithm of the concentration parameter.
%     SUBJ      – struct with field ``y`` containing a 1×16 vector of
%                 Likert responses (1…6) or NaN.
%
%   Output:
%     LL        – scalar log‑likelihood of the observed responses under
%                 the baseline model.
%
%   See also MENTALIZING_MODEL_EXP3, LEVEL0_MODEL_EXP3, LEVEL1_MODEL_EXP3,
%            LEVEL2_MODEL_EXP3.

ce = exp(parameter(1));
responses = subj.y;
w = 0.5;  % symmetric weight for the baseline
alpha = 1 + w * ce;
beta  = 1 + (1 - w) * ce;
p = diff(betainc((0:6)/6, alpha, beta));
p = max(p, realmin);
p = p / sum(p);

ll = 0;
for i = 1:numel(responses)
    y = responses(i);
    if isnan(y) || y < 1 || y > 6
        continue;
    end
    ll = ll + log(p(round(y)));
end
end