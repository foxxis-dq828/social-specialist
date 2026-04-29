function ll = baseline_model(~, subj)
%BASELINE_MODEL  Log‑likelihood for the baseline (uniform) model.
%
%   LL = BASELINE_MODEL(~, SUBJ) returns the log‑likelihood for a subject's
%   responses under a uniform model where each of the six Likert responses
%   (1…6) is equally likely.  This function ignores its first input
%   parameter (retained for compatibility with cbm_lap/cbm_hbi) and uses
%   only the response vector SUBJ.y.

responses = subj.y;
K = 6;  % number of Likert categories

% Each response contributes log(1/6) unless it is missing (NaN)
ll = sum(~isnan(responses)) * log(1/K);
end