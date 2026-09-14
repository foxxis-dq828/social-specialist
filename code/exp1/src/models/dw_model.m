function ll = dw_model(parameter, subj)
%DW_MODEL  Log-likelihood for the DW model.
%
%   LL = DW_MODEL(PARAMETER, SUBJ) computes the log-likelihood for a
%   subject's responses under the DW model. PARAMETER is a 1x3 vector:
%   PARAMETER(1) = log(ce) sets the precision,
%   PARAMETER(2) = log(n) controls mixing, and PARAMETER(3) = log(e) is an
%   elasticity term that alters the weight of agreement from a pair of
%   agents with different knowledge

ce = exp(parameter(1)); % use a exp() transform to turn the parameter into a positive value
n  = exp(parameter(2));
e  = exp(parameter(3));

responses = subj.y;
ll = 0;
for i = 1:numel(responses)
    likert = responses(i);
    p = dw_bin_probs(ce, n, e, i);
    if ~isnan(likert) && likert >= 1 && likert <= numel(p)
        ll = ll + log(max(p(likert), realmin));
    end
end
end

function p = dw_bin_probs(ce, n, e, itemIdx)
%DW_BIN_PROBS  Probability vector over 6 Likert bins for the DW model.

% weight lists; entries correspond to eight experimental conditions
denom_n  = 2 + 2 * n;
denom_ne = 2 + 2 * n * e;
w_list2 = [0.5,          1/denom_n,      1/denom_ne, ...
           1/denom_n,  (1+2*n*e)/denom_ne, (1+2*n)/denom_n, ...
           0.5,          (1+2*n)/denom_n];

w_list1 = [0.5,          (1+2*n)/denom_n, (1+2*n*e)/denom_ne, ...
           (1+2*n)/denom_n,  1/denom_ne,     1/denom_n, ...
           0.5,          1/denom_n];

% check item index
if itemIdx < 1 || itemIdx > numel(w_list1)
    error('dw_bin_probs: itemIdx out of range (1-%d)', numel(w_list1));
end

w1 = w_list1(itemIdx);
w2 = w_list2(itemIdx);

% convert DW weights into beta parameters for the predicted belief distribution
alpha = 1 + w2 * ce;
beta  = 1 + w1 * ce;

% discretize the beta distribution into six Likert bins
F = betainc((0:6)/6, alpha, beta);
p = diff(F);

% protect and normalise
p = max(p, realmin);
p = p / sum(p);
end
