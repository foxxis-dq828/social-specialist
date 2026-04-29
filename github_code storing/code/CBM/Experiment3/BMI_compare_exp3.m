%BMI_COMPARE_EXP3  Cluster participants into strategies and compare metrics across strategy clusters (Experiment 3).
%
%   This script loads subject‑level metrics and the HBI responsibility
%   matrix for Experiment 3, assigns each subject to the model with maximum
%   responsibility (marking those with responsibility below confThr as
%   "Uncertain"), optionally excludes uncertain assignments, and performs
%   one‑way ANOVA and non‑parametric tests on several metrics.  It also
%   computes group descriptives and writes the results to CSV files.
%
%   See also BMI_COMPARE_EXP1, BMI_COMPARE_EXP2.

clear; clc;

% -------- Paths --------
data_file = '/Users/foxxis/Documents/social induction/experiment data/exp3/final_data.xlsx';
hbi_file  = 'hbi_compare_exp3.mat';

% -------- Options --------
confThr       = 0.6;
excludeUncert = true;
% List of metric variables expected in the data file
needVars      = {'div','nodiv1','nodiv2','t','div_t','nodiv_t'};

% -------- Load metrics table --------
T = readtable(data_file);
missingVars = setdiff(needVars, T.Properties.VariableNames);
assert(isempty(missingVars), 'CSV missing columns: %s', strjoin(missingVars, ', '));

% Ensure numeric (robust to string/cell import)
for i = 1:numel(needVars)
    vn = needVars{i};
    if ~isnumeric(T.(vn))
        T.(vn) = str2double(string(T.(vn)));
    end
end

% -------- Load HBI responsibility --------
S = load(hbi_file, 'cbm');
assert(isfield(S, 'cbm') && isfield(S.cbm, 'output') && isfield(S.cbm.output, 'responsibility'), ...
    'No cbm.output.responsibility in %s.', hbi_file);
cbm  = S.cbm;
resp = cbm.output.responsibility;

% Model names
if isfield(cbm, 'input') && isfield(cbm.input, 'model_names') && ~isempty(cbm.input.model_names)
    model_names = cellstr(cbm.input.model_names);
else
    model_names = strcat('M', string(1:size(resp,2)));
    model_names = cellstr(model_names);
end

% Ensure resp is Nsubj × Nmodels
if size(resp,2) ~= numel(model_names) && size(resp,1) == numel(model_names)
    resp = resp.';
end
assert(size(resp,2) == numel(model_names), 'Responsibility dimension mismatch.');

% Row alignment check
assert(height(T) == size(resp,1), 'Row count mismatch: CSV=%d, responsibility=%d.', height(T), size(resp,1));

% Normalise if needed
rs = sum(resp,2);
if max(abs(rs-1)) > 1e-3
    resp = resp ./ sum(resp,2);
end

% -------- Assign strategy --------
[respMax, idx] = max(resp, [], 2);
strategy = string(model_names(idx));
strategy(respMax < confThr) = 'Uncertain';
T.strategy = categorical(strategy(:));
T.respMax  = respMax(:);

% Drop rows with missing metrics after strategy assignment
T = rmmissing(T, 'DataVariables', needVars);

if excludeUncert
    T = T(T.strategy ~= categorical({'Uncertain'}), :);
end
T.strategy = removecats(T.strategy);

% Group counts
[counts, gnames] = groupcounts(T.strategy);
disp('--- Group counts (after filtering) ---');
disp(table(cellstr(gnames), counts, 'VariableNames', {'strategy','N'}));

T1=T;

%% Mixed Effect ANOVA
Y = T(:, {'div','nodiv1','nodiv2'});

Between = table(categorical(T.strategy), 'VariableNames', {'strategyU'});

T = [Between Y];

T.strategyU = removecats(categorical(T.strategyU));


%% 
within = table(categorical({'div';'nodiv1';'nodiv2'}), 'VariableNames', {'Condition'});

rm = fitrm(T, 'div,nodiv1,nodiv2 ~ strategyU', 'WithinDesign', within);

ranovatbl = ranova(rm, 'WithinModel', 'Condition');

betweentbl = anova(rm);

disp("=== Within-subject (Condition) + Interaction (Condition*strategyU) ===")
disp(ranovatbl)

disp("=== Between-subject main effect (strategyU) ===")
disp(betweentbl)


mc = multcompare(rm, 'Condition', ...
                 'By', 'strategyU', ...              % 按 strategyU 分层
                 'ComparisonType', 'bonferroni');    % 多重比较校正（可选）

disp(mc)

mc2 = multcompare(rm, 'strategyU', ...
                  'By', 'Condition', ...
                  'ComparisonType', 'tukey-kramer');

disp(mc2)


% -------- Stats per metric --------
T=T1;

Results = table('Size',[0 8], ...
    'VariableTypes', {'string','double','double','double','double','double','double','double'}, ...
    'VariableNames', {'metric','N','kGroups','F','p_anova','eta2','p_levene','p_kw'});

PosthocAll = table();

for k = 1:numel(needVars)
    m = needVars{k};
    y = T.(m);
    g = T.strategy;

    ok = ~isnan(y) & ~isundefined(g);
    y2 = y(ok);
    g2 = removecats(g(ok));

    [gID, gNames] = findgroups(g2);
    K = numel(gNames);
    if K < 2, continue; end

    % Levene test (median)
    med = splitapply(@median, y2, gID);
    dev = abs(y2 - med(gID));
    pLev = anova1(dev, gID, 'off');

    % ANOVA
    [pA, tblA, stats] = anova1(y2, gID, 'off');
    Fstat = tblA{2,5};
    SSb  = tblA{2,2};
    SSt  = tblA{4,2};
    eta2 = SSb / SSt;

    % Kruskal‑Wallis
    pKW = kruskalwallis(y2, gID, 'off');

    % Tukey‑Kramer posthoc
    mc = multcompare(stats, 'CType','tukey-kramer', 'Display','off');
    lab = cellstr(gNames);
    ph = table();
    ph.metric = repmat(string(m), size(mc,1), 1);
    ph.group1 = string(lab(mc(:,1)));
    ph.group2 = string(lab(mc(:,2)));
    ph.diff   = mc(:,4);
    ph.CI_low = mc(:,3);
    ph.CI_high= mc(:,5);
    ph.p      = mc(:,6);
    PosthocAll = [PosthocAll; ph]; %#ok<AGROW>

    Results = [Results; {string(m), numel(y2), K, Fstat, pA, eta2, pLev, pKW}]; %#ok<AGROW>
end

disp('=== ANOVA Summary ===');
disp(Results);
disp('=== Posthoc (Tukey‑Kramer) ===');
disp(PosthocAll);

% -------- Group descriptives --------
g = T.strategy;
[Gid, gNames] = findgroups(g);
Desc = table(cellstr(gNames), accumarray(Gid,1), 'VariableNames', {'strategy','N'});

for k = 1:numel(needVars)
    m  = needVars{k};
    y  = T.(m);
    mu = splitapply(@(x) mean(x,'omitnan'), y, Gid);
    sd = splitapply(@(x) std(x,'omitnan'),  y, Gid);
    Desc.(['mean_' m]) = mu;
    Desc.(['sd_'  m]) = sd;
end

disp('=== Group descriptives (mean, SEM) ===');
disp(Desc);

