function ll = mentalizing_model(parameter, subj)
%MENTALIZING_MODEL  Log‑likelihood for the mentalizing model.
%
%   LL = MENTALIZING_MODEL(PARAMETER, SUBJ) returns the sum of log‑likelihoods
%   for a single subject's 1×8 Likert response vector under the mentalizing
%   model.  PARAMETER should be a 1×1 vector containing the logarithm of the
%   precision parameter CE (i.e., PARAMETER(1) = log(ce)).  SUBJ is a
%   structure with field ``y`` containing a 1×8 vector of integer ratings in
%   the range 1…6 (use NaN for missing data).  The model assumes that
%   responses arise from a Beta distribution whose weight ``w`` depends on
%   the item index.  See MENTALIZING_BIN_PROBS for details on the weight list.

% convert log_ce to the positive scale
ce = exp(parameter(1));

responses = subj.y;
ll = 0;

for i = 1:numel(responses)
    likert = responses(i);
    % compute bin probabilities for the i‑th item
    p = mentalizing_bin_probs(ce, i);
    if ~isnan(likert) && likert >= 1 && likert <= numel(p)
        % accumulate log‑likelihood, protecting against log(0)
        ll = ll + log(max(p(likert), realmin));
    end
end
end

function p = mentalizing_bin_probs(ce, itemIdx)
%MENTALIZING_BIN_PROBS  Probability vector over 6 Likert bins.
%
%   P = MENTALIZING_BIN_PROBS(CE, ITEMIDX) returns a 1×6 vector of
%   probabilities for Likert response bins (1…6).  CE is the precision
%   parameter on the positive scale and ITEMIDX specifies which item of
%   eight possible conditions is being evaluated.  The weight list below
%   comes from the experimental design for experiment 1.  Modify it to
%   reflect your experimental setup.

w_list_exp1 = [0.8000000, 0.2380952, 0.0625000, 0.4571429, 0.8000000, 0.8000000, 0.8000000, 0.8000000];

% guard against index errors
if itemIdx < 1 || itemIdx > numel(w_list_exp1)
    error('mentalizing_bin_probs: itemIdx out of range (1–%d)', numel(w_list_exp1));
end

w    = w_list_exp1(itemIdx);
alpha = 1 + w * ce;
beta  = 1 + (1 - w) * ce;

% compute cumulative probabilities at bin edges 0, 1/6, …, 1
cdfVals = betacdf((0:6)/6, alpha, beta);
p = diff(cdfVals);

% normalise and protect against underflow
p = max(p, realmin);
p = p / sum(p);
end