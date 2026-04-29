function ll = level0_model(parameter, subj)
%LEVEL0_MODEL  Log‑likelihood for the level‑0 model.
%
%   LL = LEVEL0_MODEL(PARAMETER, SUBJ) computes the log‑likelihood for a
%   subject's responses under the level‑0 model.  PARAMETER is a 1×2 vector
%   where PARAMETER(1) = log(ce) controls the precision of the Beta
%   distribution and PARAMETER(2) = log(n) controls a mixing parameter that
%   influences the weight lists (see LEVEL0_BIN_PROBS).  SUBJ is a structure
%   containing a field ``y`` with a 1×8 vector of Likert ratings (integers
%   between 1 and 6 or NaN).

ce = exp(parameter(1));
n  = exp(parameter(2));

responses = subj.y;
ll = 0;

for i = 1:numel(responses)
    likert = responses(i);
    p = level0_bin_probs(ce, n, i);
    if ~isnan(likert) && likert >= 1 && likert <= numel(p)
        ll = ll + log(max(p(likert), realmin));
    end
end
end

function p = level0_bin_probs(ce, n, itemIdx)
%LEVEL0_BIN_PROBS  Probability vector over 6 Likert bins for the level‑0 model.
%
%   P = LEVEL0_BIN_PROBS(CE, N, ITEMIDX) returns a 1×6 vector of bin
%   probabilities.  CE and N are positive scalars controlling the shape of
%   the Beta distribution.  The item index ITEMIDX selects the appropriate
%   weight parameters for the condition.  See the associated publication
%   for details on the derivation of the w_list1 and w_list2 vectors.

% compute denominators once for efficiency
denom = 2 + 2 * n;

% weight lists defined for experiment 1
w_list2 = [0.5,          1/denom,      1/denom,      1/denom, ...
           (1+2*n)/denom,(1+2*n)/denom,0.5,          (1+2*n)/denom];
w_list1 = [0.5,          (1+2*n)/denom,(1+2*n)/denom,(1+2*n)/denom, ...
           1/denom,      1/denom,      0.5,          1/denom];

% check item index
if itemIdx < 1 || itemIdx > numel(w_list1)
    error('level0_bin_probs: itemIdx out of range (1–%d)', numel(w_list1));
end

w1 = w_list1(itemIdx);
w2 = w_list2(itemIdx);

alpha = 1 + w2 * ce;
beta  = 1 + w1 * ce;

F = betainc((0:6)/6, alpha, beta);
p = diff(F);

% normalise and protect
p = max(p, realmin);
p = p / sum(p);
end