function plot_choice_distribution_beta_exp2(data, cbm, bestModelPerSubj)
%PLOT_CHOICE_DISTRIBUTION_BETA_EXP2  Plot choice distributions and beta predictions (Experiment 2).
%
%   PLOT_CHOICE_DISTRIBUTION_BETA_EXP2(DATA, CBM, BESTMODELPERSUBJ) creates a
%   series of bar charts showing the distribution of actual choices (1…6)
%   made by participants for each trial and strategy, overlaid with the
%   predicted distribution derived from the mean parameters of that
%   strategy (Experiment 2).  Only strategies with at least six “certain” participants
%   (max responsibility ≥ 0.5) are included.  The function lays out
%   nGroups × 8 subplots, where nGroups is the number of qualifying
%   strategies, and each subplot corresponds to one trial.  Within each
%   subplot, bars represent the count of responses in each Likert category,
%   and a red line shows the predicted beta distribution scaled to the
%   number of observations.

% Configuration
confThr   = 0.50;
minGroupN = 6;

% Strategy names
if isfield(cbm, 'input') && isfield(cbm.input, 'model_names') && ~isempty(cbm.input.model_names)
    strategyNames = cellstr(cbm.input.model_names);
else
    strategyNames = {'Mentalizing','Level0','DiversityHeuristic','ColourModel','Baseline'};
end
nModels = numel(strategyNames);

% Number of subjects and trials
N = numel(data);
T = numel(data{1}.y);

% Collect actual responses into an N×T matrix
choicesMat = nan(N, T);
for s = 1:N
    choicesMat(s,:) = data{s}.y;
end

% Responsibilities and assignments
resp = cbm.output.responsibility;
[respMax, stratIdx] = max(resp, [], 2);
validIdx = respMax >= confThr;

% Group sizes for valid subjects only
groupSizes = accumarray(stratIdx(validIdx), 1, [nModels, 1]);

% Determine which strategies have enough participants
groupsToPlot = find(groupSizes >= minGroupN);
if isempty(groupsToPlot)
    warning('No strategy groups with at least %d subjects to plot.', minGroupN);
    return;
end

% Create a figure with subplots: nGroups × T
nGroups = numel(groupsToPlot);
figure;
plotIdx = 1;

for gi = 1:nGroups
    modelIdx = groupsToPlot(gi);
    groupMask = (stratIdx == modelIdx) & validIdx;
    grpIdx    = find(groupMask);
    nGrp      = numel(grpIdx);
    % Mean parameter vector across group (in log space)
    paramMat  = cbm.output.parameters{modelIdx};
    if isempty(paramMat)
        meanParams = [];
    else
        meanParams = nanmean(paramMat(grpIdx, :), 1);
    end
    % Loop through trials
    for t = 1:T
        subplot(nGroups, T, plotIdx);
        plotIdx = plotIdx + 1;
        % Actual responses for this trial and group
        y = choicesMat(grpIdx, t);
        y = y(~isnan(y));
        % Histogram counts for categories 1..6
        counts = histcounts(y, 0.5:1:6.5);
        bar(1:6, counts, 'FaceColor',[0.6 0.6 0.9], 'EdgeColor','none');
        hold on;
        % Predicted probabilities for this trial using mean parameters
        pPred = compute_predicted_bin_probs_exp2(modelIdx, meanParams, t);
        % Scale to expected counts (same total as actual)
        predCounts = pPred * numel(y);
        plot(1:6, predCounts, 'r-', 'LineWidth',2);
        hold off;
        % Aesthetics
        ylim([0, max([counts, predCounts]) * 1.1 + 1]);
        xlim([0.5,6.5]);
        xticks(1:6);
        if gi == 1
            title(sprintf('Trial %d', t));
        end
        if t == 1
            ylabel(sprintf('%s\nCount', strategyNames{modelIdx}));
        end
        if gi == nGroups && t == ceil(T/2)
            xlabel('Rating');
        end
    end
end

sgtitle('Choice Distributions and Predicted Beta Distributions by Strategy and Trial (Experiment 2)');
end
