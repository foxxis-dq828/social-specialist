function ll = kw_model(parameter, subj)
%KW_MODEL  Log-likelihood for the KW / level-2 model.
%
%   LL = kw_model(parameter, subj) computes the log-likelihood for a
%   subject's responses under the Experiment 3 KW / level-2 model.
%   PARAMETER is a 1x2 vector where PARAMETER(1) = log(ce) and
%   PARAMETER(2) = log(n).

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

% colour gains for ear and foot conditions in experiment 3
eargain  = 1.675676;
footgain = 1.027027;

% weight list of votes defined for experiment 3
w_list = [ ...
    1 / (2 + (eargain + footgain) * n), ...
    (1 + (eargain + footgain) * n) / (2 + (eargain + footgain) * n), ...
    (1 + footgain * n) / (2 + (eargain + footgain) * n), ...
    (1 + eargain * n) / (2 + (eargain + footgain) * n), ...
    1 / (2 + eargain * 2 * n), ...
    (1 + eargain * 2 * n) / (2 + eargain * 2 * n), ...
    1 / (2 + footgain * 2 * n), ...
    (1 + footgain * 2 * n) / (2 + footgain * 2 * n), ...
    1 / (2 + (eargain + footgain) * n), ...
    (1 + (eargain + footgain) * n) / (2 + (eargain + footgain) * n), ...
    (1 + eargain * n) / (2 + (eargain + footgain) * n), ...
    (1 + footgain * n) / (2 + (eargain + footgain) * n), ...
    1 / (2 + eargain * 2 * n), ...
    (1 + eargain * 2 * n) / (2 + eargain * 2 * n), ...
    1 / (2 + footgain * 2 * n), ...
    (1 + footgain * 2 * n) / (2 + footgain * 2 * n) ...
    ];

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
