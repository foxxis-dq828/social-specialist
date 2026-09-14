function ll = baseline_model(parameter, subj)
%BASELINE_MODEL  Log-likelihood for the baseline (uniform) model.
%
% LL =  baseline_model(~, subj) returns the log-likelihood for a subject's
% responses under a uniform model where each of the six Likert responses
% (1...6) is equally likely.  
% The baseline model had no effective free parameters; a dummy parameter was 
% included only to satisfy the CBM optimization interface and 
% did not enter the likelihood.

responses = subj.y;
K = 6;  % number of Likert categories

% Each response contributes log(1/6) unless it is missing (NaN)
ll = sum(~isnan(responses)) * log(1/K);
end
