function plot_overall_means_exp3(data, cbm, bestModelPerSubj)
%PLOT_OVERALL_MEANS  Plot overall mean ratings and predicted means.
%
%   PLOT_OVERALL_MEANS(DATA, CBM, BESTMODELPERSUBJ) computes the mean
%   Likert rating for each of the eight trials across all participants and
%   displays them with a shaded error band representing the standard error
%   of the mean.  It also computes the predicted mean rating for each
%   participant on each trial based on their best‑fit model and parameter
%   values (provided by CBM) and plots the average predicted mean as a
%   dashed line.  No error bar is drawn for the predicted mean.  DATA is
%   a cell array of subject structures with field ``y`` (1×8 vector), CBM
%   is the HBI result structure, and BESTMODELPERSUBJ is an N×1 vector
%   containing the index of the best model for each subject.

% Number of subjects and trials
N = numel(data);
T = numel(data{1}.y);

% Collect actual responses into an N×T matrix
choicesMat = nan(N, T);
for s = 1:N
    choicesMat(s,:) = data{s}.y;
end

% Compute mean and standard error across subjects for each trial
meanChoices = nanmean(choicesMat, 1);
seChoices   = nanstd(choicesMat, 0, 1) ./ sqrt(max(sum(~isnan(choicesMat),1),1));

% Compute predicted mean ratings per subject and trial
predictedMat = nan(N, T);
for s = 1:N
    modelIdx = bestModelPerSubj(s);
    params   = cbm.output.parameters{modelIdx}(s, :);
    for t = 1:T
        predictedMat(s, t) = compute_predicted_mean_exp3(modelIdx, params, t);
    end
end
predictedMean = nanmean(predictedMat, 1);

% Plot
figure;
hold on;
x = 1:T;

% Shaded error region for actual mean ± SEM
upper = meanChoices + seChoices;
lower = meanChoices - seChoices;
fill([x, fliplr(x)], [upper, fliplr(lower)], [0.8 0.8 1], ...
    'EdgeColor','none', 'FaceAlpha',0.4);

% Plot actual mean (solid)
plot(x, meanChoices, 'b-', 'LineWidth',2);

% Plot predicted mean (dashed)
plot(x, predictedMean, 'r--', 'LineWidth',2);

% Labels and legend
xlabel('Trial');
ylabel('Mean Likert rating');
legend({'±SEM','Actual mean','Predicted mean'}, 'Location','best');
title('Overall Mean Ratings vs Predicted Mean');

hold off;
end