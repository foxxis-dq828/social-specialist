function m = compute_predicted_mean(modelIdx, params, trialIdx)
%COMPUTE_PREDICTED_MEAN  Expected Likert rating under a behavioural model.
%
%   M = COMPUTE_PREDICTED_MEAN(MODELIDX, PARAMS, TRIALIDX) returns the
%   expected rating on a 6‑point Likert scale (values 1…6) for a given
%   behavioural model and parameter vector.  MODELIDX specifies which
%   model is used (1=mentalizing, 2=level‑0, 3=diversity heuristic, 4=colour
%   model, 5=baseline).  PARAMS is the subject‑specific parameter vector
%   (in log space) and TRIALIDX is the trial index (1…8).  This function
%   calls COMPUTE_PREDICTED_BIN_PROBS to obtain the probability mass
%   function over the six categories and then computes the expectation.

p = compute_predicted_bin_probs(modelIdx, params, trialIdx);
% Likert categories 1 through 6
categories = 1:6;
m = sum(categories .* p);
end