%% Plot custom posterior predictions for the five Exp1 models
%
% Change only the parameter values in this first block, then run the script.
% Parameters are entered on the positive scale; the script converts them to
% log space internally to match src/models/*_model.m.

clear; clc;

projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(projectRoot, 'src', 'models'));

% Set these values yourself.
paramValues.baseline    = [];
paramValues.mentalizing = 5;          % [ce]
paramValues.vt          = [5, 1];     % [ce, n]
paramValues.dw          = [5, 1, 1];  % [ce, n, e]
paramValues.kw          = [5, 1];     % [ce, n]

modelNames = {'baseline', 'mentalizing', 'vt', 'dw', 'kw'};
modelFiles = {'baseline_model.m', 'mentalizing_model.m', 'vt_model.m', ...
              'dw_model.m', 'kw_model.m'};
barModelNames = {'mentalizing', 'vt', 'dw', 'kw'};
barModelLabels = {'M', 'VT', 'DW', 'KW'};
barColors = [0.96 0.88 0.46;   % M
             0.55 0.77 0.65;   % VT
             0.82 0.72 0.84;   % DW
             0.64 0.76 0.92];  % KW
trialsToPlot = 1:8;

for m = 1:numel(modelFiles)
    assert(exist(modelFiles{m}, 'file') == 2, ...
        'Cannot find %s. Check that src/models is on the MATLAB path.', modelFiles{m});
end

nTrials = 8;
nBins = 6;
posterior = struct();
posteriorWeight = struct();

for m = 1:numel(modelNames)
    name = modelNames{m};
    params = paramValues.(name);
    if isempty(params)
        logParams = [];
    else
        validateattributes(params, {'numeric'}, {'row', 'positive', 'finite'}, ...
            mfilename, sprintf('paramValues.%s', name));
        logParams = log(params);
    end

    posterior.(name) = nan(nTrials, nBins);
    posteriorWeight.(name) = nan(nTrials, 1);
    for trialIdx = 1:nTrials
        [posterior.(name)(trialIdx, :), posteriorWeight.(name)(trialIdx)] = ...
            model_trial_posterior(name, logParams, trialIdx);
    end
end

%% Bar plot: posterior weight for each trial
%
% For VT and DW, this is exactly w_list2.
% For mentalizing, this is the model's w for each trial.
% For KW, this is the weight used on the alpha side of the beta posterior.
posteriorBarMat = nan(numel(trialsToPlot), numel(barModelNames));
for m = 1:numel(barModelNames)
    posteriorBarMat(:, m) = posteriorWeight.(barModelNames{m})(trialsToPlot);
end
% Show the complementary probability for Trials 1 and 2 only.
complementRows = ismember(trialsToPlot, [1 2]);
posteriorBarMat(complementRows, :) = 1 - posteriorBarMat(complementRows, :);

figure('Color', [0.95 0.98 0.94], 'Name', 'Model Prediction');
tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for panelIdx = 1:numel(trialsToPlot)
    trialIdx = trialsToPlot(panelIdx);
    nexttile;
    b = bar(posteriorBarMat(panelIdx, :), 0.78, 'FaceColor', 'flat', ...
        'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 1.1);
    b.CData = barColors;

    ylim([0 1]);
    xlim([0.4 numel(barModelLabels) + 0.6]);
    xticks(1:numel(barModelLabels));
    xticklabels(barModelLabels);
    yticks(0:0.2:1);
    if ismember(trialIdx, [1 2])
        ylabel('1 - P(Type II | s1,s2)');
    else
        ylabel('P(Type II | s1,s2)');
    end
    title(sprintf('Trial %d', trialIdx), 'FontSize', 10, 'FontWeight', 'normal');
    box on;
    set(gca, 'Color', 'w', 'FontSize', 11, 'LineWidth', 1.1, ...
        'TickDir', 'out', 'XGrid', 'off', 'YGrid', 'off');
end
sgtitle('Model Prediction', 'FontWeight', 'bold', 'FontSize', 12);

%% Optional heatmap: P(response category | model, parameters, trial)
figure('Color', 'w', 'Name', 'Custom model posteriors');
tiledlayout(numel(modelNames), 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for m = 1:numel(modelNames)
    name = modelNames{m};
    nexttile;
    imagesc(posterior.(name));
    caxis([0 1]);
    colormap(parula);
    colorbar;
    xlabel('Likert response');
    ylabel('Trial');
    title(sprintf('%s model', name), 'Interpreter', 'none');
    xticks(1:nBins);
    yticks(1:nTrials);

    for trialIdx = 1:nTrials
        for binIdx = 1:nBins
            text(binIdx, trialIdx, sprintf('%.2f', posterior.(name)(trialIdx, binIdx)), ...
                'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', 'w');
        end
    end
end

%% Optional line plot: P(response category | model, parameters, trial)
figure('Color', 'w', 'Name', 'Custom model posteriors by trial');
tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

for trialIdx = 1:nTrials
    nexttile;
    hold on;
    for m = 1:numel(modelNames)
        name = modelNames{m};
        plot(1:nBins, posterior.(name)(trialIdx, :), '-o', 'LineWidth', 1.5);
    end
    hold off;
    ylim([0 1]);
    xlim([1 nBins]);
    xticks(1:nBins);
    xlabel('Likert response');
    ylabel('Posterior probability');
    title(sprintf('Trial %d', trialIdx));
    grid on;
    if trialIdx == 1
        legend(modelNames, 'Interpreter', 'none', 'Location', 'best');
    end
end

%% Local helper functions
function [p, posteriorWeight] = model_trial_posterior(modelName, logParams, trialIdx)
switch modelName
    case 'baseline'
        p = ones(1, 6) / 6;
        posteriorWeight = 0.5;

    case 'mentalizing'
        ce = exp(logParams(1));
        w_list = [0.8000000, 0.2380952, 0.0625000, 0.4571429, ...
                  0.7999999, 0.8000000, 0.7999999, 0.8000000];
        w = w_list(trialIdx);
        alpha = 1 + w * ce;
        beta = 1 + (1 - w) * ce;
        p = beta_bin_posterior(alpha, beta);
        posteriorWeight = w;

    case 'vt'
        ce = exp(logParams(1));
        n = exp(logParams(2));
        denom = 2 + 2 * n;
        w_list2 = [0.5, 1/denom, 1/denom, 1/denom, ...
                   (1 + 2*n)/denom, (1 + 2*n)/denom, 0.5, (1 + 2*n)/denom];
        w_list1 = [0.5, (1 + 2*n)/denom, (1 + 2*n)/denom, (1 + 2*n)/denom, ...
                   1/denom, 1/denom, 0.5, 1/denom];
        posteriorWeight = w_list2(trialIdx);
        alpha = 1 + posteriorWeight * ce;
        beta = 1 + w_list1(trialIdx) * ce;
        p = beta_bin_posterior(alpha, beta);

    case 'dw'
        ce = exp(logParams(1));
        n = exp(logParams(2));
        e = exp(logParams(3));
        denom_n = 2 + 2 * n;
        denom_ne = 2 + 2 * n * e;
        w_list2 = [0.5, 1/denom_n, 1/denom_ne, ...
                   1/denom_n, (1 + 2*n*e)/denom_ne, (1 + 2*n)/denom_n, ...
                   0.5, (1 + 2*n)/denom_n];
        w_list1 = [0.5, (1 + 2*n)/denom_n, (1 + 2*n*e)/denom_ne, ...
                   (1 + 2*n)/denom_n, 1/denom_ne, 1/denom_n, ...
                   0.5, 1/denom_n];
        posteriorWeight = w_list2(trialIdx);
        alpha = 1 + posteriorWeight * ce;
        beta = 1 + w_list1(trialIdx) * ce;
        p = beta_bin_posterior(alpha, beta);

    case 'kw'
        ce = exp(logParams(1));
        n = exp(logParams(2));
        d_blue = [1.870504, 1.108434, 1.108434, 1.108434];
        d_yellow = [3.085106, 1.870504, 1.108434, 1.870504];
        w_list = zeros(1, 8);
        w_list(1) = (1 + d_yellow(1)*n) / (2 + (d_blue(1) + d_yellow(1))*n);
        w_list(2) = 1 / (2 + d_blue(1)*2*n);
        w_list(3) = 1 / (2 + (d_blue(2) + d_yellow(2))*n);
        w_list(4) = 1 / (2 + d_blue(2)*2*n);
        w_list(5) = (1 + (d_blue(3) + d_yellow(3))*n) / ...
                    (2 + (d_blue(3) + d_yellow(3))*n);
        w_list(6) = (1 + d_blue(3)*2*n) / (2 + d_blue(3)*2*n);
        w_list(7) = (1 + d_blue(4)*n) / (2 + (d_blue(4) + d_yellow(4))*n);
        w_list(8) = (1 + d_blue(4)*2*n) / (2 + d_blue(4)*2*n);
        w1 = w_list(trialIdx);
        posteriorWeight = w1;
        w2 = 1 - w1;
        alpha = 1 + w1 * ce;
        beta = 1 + w2 * ce;
        p = beta_bin_posterior(alpha, beta);

    otherwise
        error('Unknown model name: %s', modelName);
end
end

function p = beta_bin_posterior(alpha, beta)
edges = (0:6) / 6;
F = betainc(edges, alpha, beta);
p = diff(F);
p = max(p, realmin);
p = p / sum(p);
end
