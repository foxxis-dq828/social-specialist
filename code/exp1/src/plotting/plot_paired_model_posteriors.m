%% Paired model posteriors: four pairs, four models, and within-pair contrasts
% Parameters below are on the positive scale (AFTER exp).
% Weight formulas match the Exp1 models and do not depend on ce.
clear; clc;

paramValues.baseline    = [];
paramValues.mentalizing = 5;          % [ce]
paramValues.vt          = [5, 1];     % [ce, n]
paramValues.dw          = [5, 1, 1];  % [ce, n, e]
paramValues.kw          = [5, 1];     % [ce, n]
% For fitted log parameters, use e.g. paramValues.dw = exp(fittedLogParams).
complementTrials = [1 2];
saveFigure = true;
% Same file and trial columns used by BMI_run_exp1.mlx.
behaviorFile = '/Users/foxxis/Documents/social induction/experiment data/exp1/output100.xlsx';
behaviorColumns = {'d1','d2','d3','d4','d5','d6','d7','d8'};

projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
outputDir = fullfile(projectRoot, 'figures', 'paired_posterior_matlab');
modelNames = {'mentalizing', 'vt', 'dw', 'kw'};
modelLabels = {'M', 'VT', 'DW', 'KW'};
parameterCounts = [1 2 3 2];
for k = 1:4
    validateattributes(paramValues.(modelNames{k}), {'numeric'}, ...
        {'real', 'row', 'numel', parameterCounts(k), 'positive', 'finite'});
end
validateattributes(complementTrials, {'numeric'}, ...
    {'vector', 'integer', '>=', 1, '<=', 8});

%% Internal posterior weights, before any plotting transformation
posteriorWeight = zeros(8, 4);
posteriorWeight(:, 1) = [0.8000000; 0.2380952; 0.0625000; 0.4571429; ...
                        0.7999999; 0.8000000; 0.7999999; 0.8000000];
n = paramValues.vt(2);
d = 2 + 2*n;
posteriorWeight(:, 2) = [0.5; 1/d; 1/d; 1/d; ...
    (1+2*n)/d; (1+2*n)/d; 0.5; (1+2*n)/d];
n = paramValues.dw(2);
e = paramValues.dw(3);
d = 2 + 2*n;
de = 2 + 2*n*e;
posteriorWeight(:, 3) = [0.5; 1/d; 1/de; 1/d; ...
    (1+2*n*e)/de; (1+2*n)/d; 0.5; (1+2*n)/d];
n = paramValues.kw(2);
dBlue = [1.870504, 1.108434, 1.108434, 1.108434];
dYellow = [3.085106, 1.870504, 1.108434, 1.870504];
posteriorWeight(:, 4) = [ ...
    (1+dYellow(1)*n)/(2+(dBlue(1)+dYellow(1))*n); ...
    1/(2+2*dBlue(1)*n); ...
    1/(2+(dBlue(2)+dYellow(2))*n); ...
    1/(2+2*dBlue(2)*n); ...
    (1+(dBlue(3)+dYellow(3))*n)/(2+(dBlue(3)+dYellow(3))*n); ...
    (1+2*dBlue(3)*n)/(2+2*dBlue(3)*n); ...
    (1+dBlue(4)*n)/(2+(dBlue(4)+dYellow(4))*n); ...
    (1+2*dBlue(4)*n)/(2+2*dBlue(4)*n)];
assert(all(isfinite(posteriorWeight(:))) && ...
    all(posteriorWeight(:) >= 0 & posteriorWeight(:) <= 1));
plotWeights = posteriorWeight;
plotWeights(complementTrials, :) = 1 - plotWeights(complementTrials, :);
pairDifferences = plotWeights(1:2:8, :) - plotWeights(2:2:8, :);

%% Participant contrasts on the same scale and orientation as predictions
assert(isfile(behaviorFile), 'Behavior file not found: %s', behaviorFile);
behaviorTable = readtable(behaviorFile);
assert(all(ismember(behaviorColumns, behaviorTable.Properties.VariableNames)), ...
    'Behavior file must contain the columns specified in behaviorColumns.');
behaviorRatings = table2array(behaviorTable(:, behaviorColumns));
validateattributes(behaviorRatings, {'numeric'}, {'real', '2d', 'ncols', 8, 'nonempty'});
validRatings = isfinite(behaviorRatings) & behaviorRatings >= 1 & ...
    behaviorRatings <= 6 & behaviorRatings == fix(behaviorRatings);
if any(~validRatings(:))
    warning('Excluding %d missing or invalid ratings from paired contrasts.', ...
        nnz(~validRatings));
end
behaviorProbabilities = (behaviorRatings - 0.5) / 6;
behaviorProbabilities(~validRatings) = NaN;
behaviorProbabilities(:, complementTrials) = 1 - behaviorProbabilities(:, complementTrials);
behaviorDifferences = behaviorProbabilities(:, 1:2:8) - behaviorProbabilities(:, 2:2:8);
behaviorPairCounts = sum(isfinite(behaviorDifferences), 1);
fprintf('Valid participant counts for Pairs 1-4: %d %d %d %d\n', behaviorPairCounts);

%% Publication figure
colors = [195 154 56; 67 134 117; 162 127 169; 88 127 165] / 255;
fig = figure('Color', 'w', 'Name', 'Paired model posteriors', ...
    'Units', 'pixels', 'Position', [100 50 900 1080]);
annotation(fig, 'textbox', [0.105 0.944 0.86 0.04], 'String', ...
    'Model predictions across paired trials', 'EdgeColor', 'none', ...
    'FontName', 'Arial', 'FontSize', 18, 'FontWeight', 'bold');
annotation(fig, 'textbox', [0.105 0.912 0.86 0.03], 'String', ...
    'Exp. 1 | Posterior weights at user-specified parameter settings', ...
    'EdgeColor', 'none', 'FontName', 'Arial', 'FontSize', 10, 'Color', [.4 .4 .4]);
annotation(fig, 'textbox', [0.105 0.849 0.5 0.03], 'String', ...
    'a   Posterior probability', 'EdgeColor', 'none', 'FontSize', 12, ...
    'FontName', 'Arial', 'FontWeight', 'bold');
annotation(fig, 'textbox', [0.62 0.849 0.37 0.03], 'String', ...
    'b   Within-pair contrast', 'EdgeColor', 'none', 'FontSize', 12, ...
    'FontName', 'Arial', 'FontWeight', 'bold');

% A common symmetric contrast scale also accommodates changed parameters.
allDifferences = [pairDifferences(:); behaviorDifferences(isfinite(behaviorDifferences))];
deltaLimit = max(0.8, ceil((max(abs(allDifferences)) + 0.16)*10)/10);
for pair = 1:4
    trials = [2*pair-1, 2*pair];
    bottom = 0.695 - (pair-1)*0.192;
    ax = axes(fig, 'Position', [0.105 bottom 0.40 0.115]);
    hold(ax, 'on');
    for k = 1:4
        bar(ax, k-0.17, plotWeights(trials(1), k), 0.30, ...
            'FaceColor', colors(k,:), 'EdgeColor', colors(k,:), 'LineWidth', 1.2);
        bar(ax, k+0.17, plotWeights(trials(2), k), 0.30, ...
            'FaceColor', 'w', 'EdgeColor', colors(k,:), 'LineWidth', 1.2);
    end
    styleAxes(ax, modelLabels);
    ylim(ax, [0 1]); yticks(ax, [0 0.5 1]);
    if all(ismember(trials, complementTrials))
        probabilityLabel = '1 - P(Type II | s_1,s_2)';
        deltaLabel = '\Delta(1-P)';
    elseif ~any(ismember(trials, complementTrials))
        probabilityLabel = 'P(Type II | s_1,s_2)';
        deltaLabel = '\Delta P';
    else
        probabilityLabel = 'Displayed probability';
        deltaLabel = '\Delta displayed probability';
    end
    ylabel(ax, probabilityLabel, 'Interpreter', 'tex');
    title(ax, sprintf('Pair %d   |   Trial %d vs %d', pair, trials), ...
        'FontSize', 10, 'FontWeight', 'bold');
    if pair == 1
        h1 = patch(ax, NaN, NaN, [.45 .45 .45], 'EdgeColor', [.35 .35 .35]);
        h2 = patch(ax, NaN, NaN, 'w', 'EdgeColor', [.35 .35 .35]);
        lg = legend(ax, [h1 h2], {'First trial in pair', 'Second trial in pair'}, ...
            'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 10);
        set(lg, 'Units', 'normalized', 'Position', [0.105 0.886 0.5 0.024]);
    end

    ax = axes(fig, 'Position', [0.63 bottom 0.34 0.115]);
    hold(ax, 'on');
    yline(ax, 0, '-', 'Color', [.35 .35 .35], 'LineWidth', 0.9);
    for k = 1:2
        value = pairDifferences(pair,k);
        plot(ax, [k k], [0 value], '-', 'Color', colors(k,:), 'LineWidth', 2);
        plot(ax, k, value, 'o', 'MarkerFaceColor', colors(k,:), ...
            'MarkerEdgeColor', colors(k,:), 'MarkerSize', 6);
        if abs(value) < 0.005
            label = '0.00';
        else
            label = sprintf('%+.2f', value);
        end
        if value >= -1e-6
            offset = 0.065; alignment = 'bottom';
        else
            offset = -0.065; alignment = 'top';
        end
        text(ax, k, value+offset, label, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', alignment, 'FontName', 'Arial', 'FontSize', 9);
    end
    % Spread tied values horizontally; their vertical values remain exact.
    observed = behaviorDifferences(:, pair);
    observed = observed(isfinite(observed));
    observedX = 3 * ones(size(observed));
    [~, ~, tieGroup] = unique(round(observed * 6));
    for group = 1:max([0; tieGroup])
        members = find(tieGroup == group);
        if numel(members) > 1
            spread = min(0.24, 0.018*(numel(members)-1));
            observedX(members) = 3 + linspace(-spread, spread, numel(members))';
        end
    end
    scatter(ax, observedX, observed, 12, [.25 .25 .25], 'filled', ...
        'MarkerFaceAlpha', 0.45, 'MarkerEdgeColor', 'none');
    xline(ax, 2.5, ':', 'Color', [.75 .75 .75], 'LineWidth', 0.65);
    styleAxes(ax, [modelLabels(1:2) {'Data'}]);
    ylim(ax, [-deltaLimit deltaLimit]); yticks(ax, [-0.5 0 0.5]);
    yticklabels(ax, {'-0.5', '0', '+0.5'});
    ylabel(ax, deltaLabel, 'Interpreter', 'tex');
    title(ax, sprintf('Trial %d - Trial %d', trials), 'FontWeight', 'normal', ...
        'FontSize', 10, 'Color', [.4 .4 .4]);
end
annotation(fig, 'textbox', [0.105 0.039 0.87 0.04], 'String', ...
    {'M, mentalizing; VT, vote tallying; DW and KW follow the Exp1 model definitions.', ...
     'Data: one point per participant; ratings mapped to (rating - 0.5)/6; horizontal spread separates ties.'}, ...
    'EdgeColor', 'none', 'FontName', 'Arial', 'FontSize', 8, 'Color', [.4 .4 .4]);

if saveFigure
    if ~exist(outputDir, 'dir'), mkdir(outputDir); end
    exportgraphics(fig, fullfile(outputDir, 'paired_posteriors.png'), 'Resolution', 300);
    exportgraphics(fig, fullfile(outputDir, 'paired_posteriors.pdf'), 'ContentType', 'vector');
    savefig(fig, fullfile(outputDir, 'paired_posteriors.fig'));
    save(fullfile(outputDir, 'paired_posteriors.mat'), 'paramValues', ...
        'posteriorWeight', 'plotWeights', 'pairDifferences', 'complementTrials', ...
        'behaviorFile', 'behaviorColumns', 'behaviorRatings', 'behaviorProbabilities', ...
        'behaviorDifferences', 'behaviorPairCounts');
end

function styleAxes(ax, labels)
set(ax, 'FontName', 'Arial', 'FontSize', 10, 'Color', 'w', ...
    'Box', 'off', 'LineWidth', 0.65, 'TickDir', 'out', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'XGrid', 'off', 'YGrid', 'on', 'GridColor', [.85 .85 .85], 'GridAlpha', 0.4);
xlim(ax, [0.4 numel(labels)+0.6]); xticks(ax, 1:numel(labels)); xticklabels(ax, labels);
end
