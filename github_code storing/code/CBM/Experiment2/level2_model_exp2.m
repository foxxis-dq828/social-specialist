function ll = level2_model_exp2(parameter, subj)
%LEVEL2_MODEL_EXP2  Log‑likelihood for the colour (level‑2) model (Experiment 2).
%
%   LL = LEVEL2_MODEL_EXP2(PARAMETER, SUBJ) returns the log‑likelihood of
%   a subject's responses under the colour model used in Experiment 2.
%   PARAMETER is [log_ce, log_n], representing the log of the concentration
%   parameter ce and the colour‑dependence parameter n.  The model uses
%   colour constants (d_blue and d_yellow) to construct a weight list for
%   each experimental condition.
%
%   Inputs:
%     PARAMETER – [log_ce, log_n].
%     SUBJ      – struct with field y.
%
%   Output:
%     LL – scalar log‑likelihood.
%
%   See also MENTALIZING_MODEL_EXP2, LEVEL0_MODEL_EXP2, LEVEL1_MODEL_EXP2.

ce = exp(parameter(1));
n  = exp(parameter(2));
responses = subj.y;

% Colour constants for Experiment 2 (duplicated across eight conditions).
d_blue   = [1.604278, 1.604278, 1.410604, 1.410604, 1.040060, 1.040060, 1.463840, 1.463840];
d_yellow = [1.604278, 1.604278, 1.410604, 1.410604, 1.040060, 1.040060, 1.463840, 1.463840];

ll = 0;
% Precompute weight list.
w_list = zeros(1,8);
w_list(1) = (1 + d_blue(1)*2*n) / (2 + d_blue(1)*2*n);
w_list(2) = (1 + (d_blue(2) + d_yellow(2))*n) / (2 + (d_blue(2) + d_yellow(2))*n);
w_list(3) = 1 / (2 + d_yellow(3)*2*n);
w_list(4) = (1 + d_blue(4)*n) / (2 + (d_yellow(4) + d_blue(4))*n);
w_list(5) = (1 + d_blue(5)*2*n) / (2 + d_blue(5)*2*n);
w_list(6) = (1 + (d_blue(6) + d_yellow(6))*n) / (2 + (d_blue(6) + d_yellow(6))*n);
w_list(7) = 1 / (2 + d_blue(7)*2*n);
w_list(8) = (1 + d_yellow(8)*n) / (2 + (d_yellow(8) + d_blue(8))*n);

for i = 1:numel(responses)
    y = responses(i);
    if isnan(y) || y < 1 || y > 6
        continue;
    end
    w1 = w_list(i);
    w2 = 1 - w1;
    alpha = 1 + w1 * ce;
    beta  = 1 + w2 * ce;
    p = diff(betainc((0:6)/6, alpha, beta));
    p = max(p, realmin);
    p = p / sum(p);
    ll = ll + log(p(round(y)));
end
end
