function m = compute_predicted_mean_exp2(modelIdx, params, trialIdx)
%COMPUTE_PREDICTED_MEAN_EXP2  Expected Likert rating under a behavioural model (Experiment 2).
%
%   M = COMPUTE_PREDICTED_MEAN_EXP2(MODELIDX, PARAMS, TRIALIDX) returns the
%   expected rating on the 6‑point Likert scale (1…6) for a given model,
%   parameter vector and experimental condition.  It calls
%   COMPUTE_PREDICTED_BIN_PROBS_EXP2 and computes the weighted average of
%   categories.
%
%   Inputs:
%     MODELIDX – index of the model (1=mentalizing, 2=level‑0, 3=level‑1, 4=level‑2, 5=baseline).
%     PARAMS   – vector of log parameters.
%     TRIALIDX – integer in 1…8 specifying the trial/condition.
%
%   Output:
%     M – expected Likert rating.
%
%   See also COMPUTE_PREDICTED_BIN_PROBS_EXP2.

p = compute_predicted_bin_probs_exp2(modelIdx, params, trialIdx);
categories = 1:6;
m = sum(categories .* p);
end
