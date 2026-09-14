function ll = kw_model(parameter, subj)
%KW_MODEL  Log-likelihood for the KW model.
%
%   LL = kw_model(parameter, subj) computes the log-likelihood for a
%   subject's responses under the KW model. PARAMETER is a 1x2 vector where
%   PARAMETER(1) = log(ce) and PARAMETER(2) = log(n).

ce = exp(parameter(1));  % use a exp() transform to turn the parameter into a positive value
n  = exp(parameter(2));

responses = subj.y;
ll = 0;
for i = 1:numel(responses)
    likert = responses(i);
    p = kw_bin_probs(ce, n, i);
    if ~isnan(likert) && likert >= 1 && likert <= numel(p)
        ll = ll + log(max(p(likert), realmin));
    end
end
end

function p = kw_bin_probs(ce, n, itemIdx)
%KW_BIN_PROBS  Probability vector over 6 Likert bins for the KW model.

% diagnosticity of different experts (blue agents, yellow agents) across 4
% conditions
d_blue   = [1.870504, 1.108434, 1.108434, 1.108434];
d_yellow = [3.085106, 1.870504, 1.108434, 1.870504];

% weight lists of votes defined for experiment 1
w_list = zeros(1, 8);
w_list(1) = (1 + d_yellow(1)*n) / (2 + (d_blue(1) + d_yellow(1))*n);
w_list(2) = 1 / (2 + (d_blue(1))*2*n);
w_list(3) = 1 / (2 + (d_blue(2) + d_yellow(2))*n);
w_list(4) = 1 / (2 + (d_blue(2))*2*n);
w_list(5) = (1 + (d_blue(3) + d_yellow(3))*n) / (2 + (d_blue(3) + d_yellow(3))*n);
w_list(6) = (1 + (d_blue(3))*2*n) / (2 + (d_blue(3))*2*n);
w_list(7) = (1 + d_blue(4)*n) / (2 + (d_blue(4) + d_yellow(4))*n);
w_list(8) = (1 + (d_blue(4))*2*n) / (2 + (d_blue(4))*2*n);

% check item index
if itemIdx < 1 || itemIdx > numel(w_list)
    error('kw_bin_probs: itemIdx out of range (1-%d)', numel(w_list));
end

w1 = w_list(itemIdx);
w2 = 1 - w1;

% convert knowledgeability weights into beta parameters for the predicted belief distribution
alpha = 1 + w1 * ce;
beta  = 1 + w2 * ce;

% discretize the beta distribution into six Likert bins
F = betainc((0:6)/6, alpha, beta);
p = diff(F);

% protect and normalise
p = max(p, realmin);
p = p / sum(p);
end
