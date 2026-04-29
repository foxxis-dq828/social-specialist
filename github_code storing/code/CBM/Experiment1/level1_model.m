function ll = level1_model(parameter, subj)
%LEVEL1_MODEL  Log‑likelihood for the level‑1 (diversity heuristic) model.
%
%   LL = LEVEL1_MODEL(PARAMETER, SUBJ) computes the log‑likelihood for a
%   subject's responses under the diversity heuristic (level‑1) model.
%   PARAMETER is a 1×3 vector: PARAMETER(1) = log(ce) sets the precision,
%   PARAMETER(2) = log(n) controls mixing, and PARAMETER(3) = log(e) is an
%   elasticity term that alters some of the weight entries.  SUBJ is a
%   structure with field ``y`` containing a 1×8 vector of Likert ratings.

ce = exp(parameter(1));
n  = exp(parameter(2));
e  = exp(parameter(3));

responses = subj.y;
ll = 0;
for i = 1:numel(responses)
    likert = responses(i);
    p = level1_bin_probs(ce, n, e, i);
    if ~isnan(likert) && likert >= 1 && likert <= numel(p)
        ll = ll + log(max(p(likert), realmin));
    end
end
end

function p = level1_bin_probs(ce, n, e, itemIdx)
%LEVEL1_BIN_PROBS  Probability vector over 6 Likert bins for the diversity heuristic.
%
%   P = LEVEL1_BIN_PROBS(CE, N, E, ITEMIDX) returns a 1×6 vector of bin
%   probabilities.  CE is the precision, N controls the weight lists, and E
%   introduces an elasticity to some weights.  The itemIdx chooses the
%   appropriate weight entries.  See the corresponding behavioural model for
%   details.

denom_n  = 2 + 2 * n;
denom_ne = 2 + 2 * n * e;

% weight lists; entries correspond to eight experimental conditions
w_list2 = [0.5,          1/denom_n,      1/denom_ne, ...
           1/denom_n,  (1+2*n*e)/denom_ne, (1+2*n)/denom_n, ...
           0.5,          (1+2*n)/denom_n];

w_list1 = [0.5,          (1+2*n)/denom_n, (1+2*n*e)/denom_ne, ...
           (1+2*n)/denom_n,  1/denom_ne,     1/denom_n, ...
           0.5,          1/denom_n];

if itemIdx < 1 || itemIdx > numel(w_list1)
    error('level1_bin_probs: itemIdx out of range (1–%d)', numel(w_list1));
end

w1 = w_list1(itemIdx);
w2 = w_list2(itemIdx);

alpha = 1 + w2 * ce;
beta  = 1 + w1 * ce;

F = betainc((0:6)/6, alpha, beta);
p = diff(F);

% protect and normalise
p = max(p, realmin);
p = p / sum(p);
end