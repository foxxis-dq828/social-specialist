function plot_choices_by_strategy(data, cbm, bestModelPerSubj)
%PLOT_CHOICES_BY_STRATEGY  Plot mean ratings by strategy group with predictions.
%
%   PLOT_CHOICES_BY_STRATEGY(DATA, CBM, BESTMODELPERSUBJ) groups
%   participants according to their best‑fit model (strategy), filters out
%   uncertain assignments (max responsibility < 0.5), and keeps only those
%   strategies with more than five participants.  For each remaining
%   strategy, the function computes the mean Likert rating for each of the
%   eight trials and plots it with a shaded error band.  It also computes
%   the predicted mean rating for each participant on each trial based on
%   their model parameters and plots the average predicted mean as a
%   dashed line.  Colours distinguish different strategies.

% Configuration: threshold for being considered “certain” and minimum group size
confThr   = 0.50;
minGroupN = 6;  % include groups with at least 6 subjects

% Strategy names in the order of cbm.output.responsibility columns
strategyNames = {'Mentalizing','Level0','DiversityHeuristic','ColourModel','Baseline'};
nModels = numel(strategyNames);

% Number of subjects and trials
N = numel(data);
T = numel(data{1}.y);

% Collect actual responses into an N×T matrix
choicesMat = nan(N, T);
for s = 1:N
    choicesMat(s,:) = data{s}.y;
end

% Compute responsibilities and assignments
resp = cbm.output.responsibility;
[respMax, stratIdx] = max(resp, [], 2);

% Filter: keep subjects with max responsibility ≥ confThr
validIdx = respMax >= confThr;

% Compute group sizes
groupSizes = accumarray(stratIdx(validIdx), 1, [nModels, 1]);

% Identify strategies meeting the minimum group size requirement
groupsToPlot = find(groupSizes >= minGroupN);
if isempty(groupsToPlot)
    warning('No strategy groups with at least %d subjects. Skipping plot.', minGroupN);
    return;
end

% Prepare colour palette
colors = lines(numel(groupsToPlot));

% Create figure
figure;
hold on;
legEntries = {};
hHandles   = [];

for gi = 1:numel(groupsToPlot)
    modelIdx = groupsToPlot(gi);
    groupMask = (stratIdx == modelIdx) & validIdx;
    if sum(groupMask) < minGroupN
        continue;
    end
    grpChoices = choicesMat(groupMask, :);
    % Compute actual mean and SEM across group
    meanChoices = nanmean(grpChoices, 1);
    seChoices   = nanstd(grpChoices, 0, 1) ./ sqrt(max(sum(~isnan(grpChoices),1),1));
    
    % Compute predicted mean ratings per subject and trial for this group
    grpIdx = find(groupMask);
    predictedMatGrp = nan(numel(grpIdx), T);
    for k = 1:numel(grpIdx)
        s = grpIdx(k);
        params = cbm.output.parameters{modelIdx}(s, :);
        for t = 1:T
            predictedMatGrp(k, t) = compute_predicted_mean(modelIdx, params, t);
        end
    end
    predictedMean = nanmean(predictedMatGrp, 1);
    
    % Plot actual mean with shaded error band
    x = 1:T;
    c = colors(gi, :);
    upper = meanChoices + seChoices;
    lower = meanChoices - seChoices;
    fill([x, fliplr(x)], [upper, fliplr(lower)], c, 'FaceAlpha',0.2, 'EdgeColor','none');
    hActual = plot(x, meanChoices, '-', 'Color', c, 'LineWidth',2);
    hPred   = plot(x, predictedMean, '--', 'Color', c, 'LineWidth',1.5);
    
    % Store handles and legend entries
    hHandles = [hHandles; hActual; hPred]; %#ok<AGROW>
    legEntries{end+1} = sprintf('%s actual (n=%d)', strategyNames{modelIdx}, sum(groupMask)); %#ok<AGROW>
    legEntries{end+1} = sprintf('%s predicted', strategyNames{modelIdx}); %#ok<AGROW>
end

% Labels and legend
xlabel('Trial');
ylabel('Mean Likert rating');
legend(hHandles, legEntries, 'Location','best', 'Interpreter','none');
title('Mean Ratings by Strategy with Predicted Means');

hold off;
end