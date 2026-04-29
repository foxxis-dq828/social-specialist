function plot_choice_distribution_beta_exp3(data, cbm, bestModelPerSubj)
%PLOT_CHOICE_DISTRIBUTION_BETA_EXP3  Plot empirical and predicted choice distributions (Experiment 3).
%
%   PLOT_CHOICE_DISTRIBUTION_BETA_EXP3(DATA, CBM, BESTMODELPERSUBJ) creates
%   bar charts for each combination of strategy (behavioural model) and
%   trial in Experiment 3, showing the distribution of actual Likert
%   responses overlaid with the predicted distribution derived from the
%   average parameter vector for that strategy.  Only strategies with at
%   least a minimum number of “certain” subjects (max responsibility ≥ 0.5)
%   are plotted.  The figure is organised into a grid of nGroups × 16
%   subplots, where nGroups is the number of strategies meeting the
%   inclusion criteria.

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
        if iscell(paramMat)
            % If parameters are stored in a cell array, extract and average
            tmp = cellfun(@(c) c(:)', paramMat(grpIdx), 'UniformOutput', false);
            meanParams = mean(cat(1, tmp{:}), 1);
        else
            meanParams = nanmean(paramMat(grpIdx, :), 1);
        end
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
        pPred = compute_predicted_bin_probs_exp3(modelIdx, meanParams, t);
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

sgtitle('Choice Distributions and Predicted Beta Distributions by Strategy and Trial (Experiment 3)');
end