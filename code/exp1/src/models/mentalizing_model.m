function ll = mentalizing_model(parameter, subj)
%MENTALIZING_MODEL  Log‑likelihood for the mentalizing model.
%
%   LL = MENTALIZING_MODEL(PARAMETER, SUBJ) returns the sum of log‑likelihoods
%   for each subject's 8 responses under the mentalizing model.  

% convert log_ce to the positive scale
ce = exp(parameter(1));

responses = subj.y;
ll = 0;

for i = 1:numel(responses)
    likert = responses(i);
    % compute bin probabilities for the i‑th item
    p = mentalizing_bin_probs(ce, i);
    % accumulate log‑likelihood, protecting against log(0)
    ll = ll + log(max(p(likert), realmin));
end
end

function p = mentalizing_bin_probs(ce, itemIdx)
%MENTALIZING_BIN_PROBS  Probability vector over 6 Likert bins.
%
%   P = mentalizing_bin_probs(ce, itemIdx) returns a 1×6 vector of
%   probabilities for Likert response bins (1…6).  ce is the precision
%   parameter on the positive scale and itemIdx specifies which condition is being evaluated.  
%   The weight list below comes from the experimental design for experiment 1. 

w_list_exp1 = [0.80000000, 0.23809524, 0.07142857, 0.45714286, 0.80000000, 0.80000000, 0.80000000, 0.80000000];

% calculate the belief distribution (beta distribution) representing
% participants' beliefs toward the identity of the jar (I/II)
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