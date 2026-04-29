function ll = level2_model(parameter, subj)
%LEVEL2_MODEL  Log‑likelihood for the level‑2/colour model.
%
%   LL = LEVEL2_MODEL(PARAMETER, SUBJ) computes the log‑likelihood for a
%   subject's responses under the level‑2 (colour) model.  PARAMETER is a
%   1×2 vector where PARAMETER(1) = log(ce) and PARAMETER(2) = log(n).
%   SUBJ is a structure with field ``y`` containing a 1×8 vector of Likert
%   ratings (1…6).  This model incorporates colour‑dependent weights
%   specified by ``d_blue`` and ``d_yellow`` vectors; see LEVEL2_BIN_PROBS.

ce = exp(parameter(1));
n  = exp(parameter(2));

responses = subj.y;
ll = 0;
for i = 1:numel(responses)
    likert = responses(i);
    p = level2_bin_probs(ce, n, i);
    if ~isnan(likert) && likert >= 1 && likert <= numel(p)
        ll = ll + log(max(p(likert), realmin));
    end
end
end

function p = level2_bin_probs(ce, n, itemIdx)
%LEVEL2_BIN_PROBS  Probability vector over 6 Likert bins for the level‑2/colour model.
%
%   P = LEVEL2_BIN_PROBS(CE, N, ITEMIDX) returns a 1×6 vector of bin
%   probabilities for a given experimental condition ITEMIDX (1…8).  CE and
%   N are positive scalars.  The weight ``w1`` is derived from the
%   colour‑specific constants ``d_blue`` and ``d_yellow`` defined below.

% colour constants for four item pairs
d_blue   = [1.870504, 1.108434, 1.108434, 1.108434];
d_yellow = [3.085106, 1.870504, 1.108434, 1.870504];

w_list = zeros(1, 8);
w_list(1) = (1 + d_yellow(1)*n) / (2 + (d_blue(1) + d_yellow(1))*n);
w_list(2) = 1 / (2 + (d_blue(1))*2*n);
w_list(3) = 1 / (2 + (d_blue(2) + d_yellow(2))*n);
w_list(4) = 1 / (2 + (d_blue(2))*2*n);
w_list(5) = (1 + (d_blue(3) + d_yellow(3))*n) / (2 + (d_blue(3) + d_yellow(3))*n);
w_list(6) = (1 + (d_blue(3))*2*n) / (2 + (d_blue(3))*2*n);
w_list(7) = (1 + d_blue(4)*n) / (2 + (d_blue(4) + d_yellow(4))*n);
w_list(8) = (1 + (d_blue(4))*2*n) / (2 + (d_blue(4))*2*n);

if itemIdx < 1 || itemIdx > numel(w_list)
    error('level2_bin_probs: itemIdx out of range (1–%d)', numel(w_list));
end

w1 = w_list(itemIdx);
w2 = 1 - w1;

alpha = 1 + w1 * ce;
beta  = 1 + w2 * ce;

F = betainc((0:6)/6, alpha, beta);
p = diff(F);

% protect and normalise
p = max(p, realmin);
p = p / sum(p);
end