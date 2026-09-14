function ll = vt_model(parameter, subj)
%VT_MODEL  Log-likelihood for the VT model.
%
%  LL = vt_model(parameter, subj) computes the log-likelihood for a
%  subject's responses under the VT model. PARAMETER is a 1x2 vector
%  where PARAMETER(1) = log(ce) controls the precision of the beta
%  distribution and PARAMETER(2) = log(n) controls the weight of each vote. 

ce = exp(parameter(1)); % use a exp() transform to turn the parameter into a positive value
n  = exp(parameter(2));

responses = subj.y;
ll = 0;

for i = 1:numel(responses)
    likert = responses(i);
    p = vt_bin_probs(ce, n, i);
    if ~isnan(likert) && likert >= 1 && likert <= numel(p)
        ll = ll + log(max(p(likert), realmin));
    end
end
end

function p = vt_bin_probs(ce, n, itemIdx)
%VT_BIN_PROBS  Probability vector over 6 Likert bins for the VT model


% weight lists defined for experiment 2
denom = 2 + 2 * n;
w_list2 = [(1+2*n)/denom, (1+2*n)/denom, 1/denom, 0.5, ...
           (1+2*n)/denom, (1+2*n)/denom, 1/denom, 0.5];
w_list1 = [1/denom,       1/denom,       (1+2*n)/denom, 0.5, ...
           1/denom,       1/denom,       (1+2*n)/denom, 0.5];

% check item index
if itemIdx < 1 || itemIdx > numel(w_list1)
    error('vt_bin_probs: itemIdx out of range (1-%d)', numel(w_list1));
end

w1 = w_list1(itemIdx);
w2 = w_list2(itemIdx);

% convert VT weights into beta parameters for the predicted belief distribution
alpha = 1 + w2 * ce;
beta  = 1 + w1 * ce;

% discretize the beta distribution into six Likert bins
F = betainc((0:6)/6, alpha, beta);
p = diff(F);

% normalise and protect
p = max(p, realmin);
p = p / sum(p);
end
