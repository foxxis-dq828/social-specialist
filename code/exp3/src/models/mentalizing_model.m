function ll = mentalizing_model(parameter, subj)
%MENTALIZING_MODEL  Log-likelihood for the mentalizing model.
%
%   LL = MENTALIZING_MODEL(PARAMETER, SUBJ) returns the sum of
%   log-likelihoods for each subject's 16 responses under the Experiment 3
%   mentalizing model.

% convert log_ce to the positive scale
ce = exp(parameter(1));

responses = subj.y;
ll = 0;

for i = 1:numel(responses)
    likert = responses(i);
    % compute bin probabilities for the i-th item
    p = mentalizing_bin_probs(ce, i);
    % accumulate log-likelihood, protecting against log(0)
    if ~isnan(likert) && likert >= 1 && likert <= numel(p)
        ll = ll + log(max(p(likert), realmin));
    end
end
end

function p = mentalizing_bin_probs(ce, itemIdx)
%MENTALIZING_BIN_PROBS  Probability vector over 6 Likert bins.
%
%   P = mentalizing_bin_probs(ce, itemIdx) returns a 1x6 vector of
%   probabilities for Likert response bins (1..6). ce is the precision
%   parameter on the positive scale and itemIdx specifies which condition is
%   being evaluated. The weight list below comes from the experimental
%   design for experiment 3.

w_list_exp3 = [0.0000000, 0.0000000, 1.0000000, 1.0000000, ...
               0.2727273, 0.7777778, 0.4666667, 0.6000000, ...
               0.0000000, 0.0000000, 1.0000000, 1.0000000, ...
               0.2727273, 0.7777778, 0.4666667, 0.6000000];
% check item index
if itemIdx < 1 || itemIdx > numel(w_list_exp3)
    error('mentalizing_bin_probs: itemIdx out of range (1-%d)', numel(w_list_exp3));
end

% calculate the belief distribution (beta distribution) representing
% participants' beliefs toward the identity of the jar (I/II)
w    = w_list_exp3(itemIdx);
alpha = 1 + w * ce;
beta  = 1 + (1 - w) * ce;

% compute cumulative probabilities at bin edges 0, 1/6, ..., 1
cdfVals = betainc((0:6)/6, alpha, beta);
p = diff(cdfVals);

% normalise and protect against underflow
p = max(p, realmin);
p = p / sum(p);
end
