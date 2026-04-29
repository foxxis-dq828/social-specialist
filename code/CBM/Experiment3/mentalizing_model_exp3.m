function ll = mentalizing_model_exp3(parameter, subj)
%MENTALIZING_MODEL_EXP3  Log‑likelihood for the mentalizing model (Experiment 3).
%
%   LL = MENTALIZING_MODEL_EXP3(PARAMETER, SUBJ) returns the log‑likelihood
%   of a subject's Likert responses under the mentalizing model used in
%   Experiment 3.  PARAMETER is a 1×1 vector containing the log of the
%   concentration parameter (log_ce).  The model assumes that ratings in
%   each of the 16 conditions are generated from a Beta distribution with
%   weight parameters drawn from a fixed list derived from the original
%   experiment.  These weights determine the mean of the Beta and are
%   multiplied by the concentration parameter ce = exp(log_ce).
%
%   Inputs:
%     PARAMETER – [log_ce], the logarithm of the concentration parameter.
%     SUBJ      – struct with field ``y`` containing a 1×16 vector of Likert
%                 responses (1…6) or NaN for missing values.
%
%   Output:
%     LL        – scalar log‑likelihood of the observed responses under
%                 the mentalizing model.
%
%   See also LEVEL0_MODEL_EXP3, LEVEL1_MODEL_EXP3, LEVEL2_MODEL_EXP3,
%            BASELINE_MODEL_EXP3.

% Exponentiate the log concentration parameter.  Clamp to avoid under/overflow.
ce = exp(parameter(1));

% Extract responses and weight list for the 16 experimental conditions.  The
% weights correspond to the subject's inferred probability of choosing the
% higher end of the scale in each condition.
responses = subj.y;
w_list = [0.0000000, 0.0000000, 1.0000000, 1.0000000, ...
          0.2735507, 0.7765487, 0.4667553, 0.5992063, ...
          0.0000000, 0.0000000, 1.0000000, 1.0000000, ...
          0.2735507, 0.7765487, 0.4667553, 0.5992063];

% Initialise log‑likelihood.
ll = 0;

for i = 1:numel(responses)
    y = responses(i);
    % Skip missing or out‑of‑range responses
    if isnan(y) || y < 1 || y > 6
        continue;
    end
    w = w_list(i);
    alpha = 1 + w * ce;
    beta  = 1 + (1 - w) * ce;
    % Compute probabilities for the six Likert bins using the incomplete beta function.
    p = diff(betainc((0:6)/6, alpha, beta));
    % Protect against zero probabilities and normalise.
    p = max(p, realmin);
    p = p / sum(p);
    % Accumulate log‑likelihood.  Use round to index into p safely.
    ll = ll + log(p(round(y)));
end
end