function ll = level2_model_exp3(parameter, subj)
%LEVEL2_MODEL_EXP3  Log‑likelihood for the colour model (Experiment 3).
%
%   LL = LEVEL2_MODEL_EXP3(PARAMETER, SUBJ) computes the log‑likelihood of
%   a subject's responses under the colour model used in Experiment 3.
%   The model has two free parameters: the log concentration (log_ce) and
%   the log colour weight (log_n).  Colour gains for ear and foot
%   conditions are defined by constants ``eargain`` and ``footgain``.
%   These gains determine a weight list for each of the 16 conditions.
%
%   Inputs:
%     PARAMETER – [log_ce, log_n].
%     SUBJ      – struct with field ``y`` containing a 1×16 vector of
%                 Likert responses (1…6) or NaN.
%
%   Output:
%     LL        – scalar log‑likelihood of the observed responses.
%
%   See also MENTALIZING_MODEL_EXP3, LEVEL0_MODEL_EXP3, LEVEL1_MODEL_EXP3,
%            BASELINE_MODEL_EXP3.

% Exponentiate parameters
ce = exp(parameter(1));
n  = exp(parameter(2));

responses = subj.y;

% Colour gains (constants) for ear and foot.  These values come from the
% original R code and are repeated across conditions as appropriate.
eargain  = 1.675676;
footgain = 1.027027;

% Construct the weight list for 16 conditions.  Each entry represents the
% probability of weighting the 'blue' option when forming the Beta mean.
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

% Initialise log‑likelihood
ll = 0;

for i = 1:numel(responses)
    y = responses(i);
    if isnan(y) || y < 1 || y > 6
        continue;
    end
    w1 = w_list(i);
    w2 = 1 - w1;
    alpha = 1 + w1 * ce;
    beta  = 1 + w2 * ce;
    p = diff(betainc((0:6)/6, alpha, beta));
    p = max(p, realmin);
    p = p / sum(p);
    ll = ll + log(p(round(y)));
end
end