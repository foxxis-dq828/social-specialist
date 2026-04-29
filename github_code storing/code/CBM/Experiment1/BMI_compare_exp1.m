%% BMI_compare_exp1.m
% Cluster participants into strategies using HBI responsibility,
% then compare av_pdiv, av_nodiv, time across strategy clusters (ANOVA).

clear; clc;

xlsx_file = '/Users/foxxis/Documents/social induction/experiment data/exp1/output100.xlsx';
hbi_file  = 'hbi_compare.mat';   % must contain variable cbm

confThr       = 0.6;   % respMax < confThr -> Uncertain

%% Load summary metrics 
Sum = readtable(xlsx_file);

needVars = {'av_pdiv','av_nodiv','time','t','median_t','div_t','nodiv_t'};
missingVars = setdiff(needVars, Sum.Properties.VariableNames);
assert(isempty(missingVars), 'Missing columns in Excel: %s', strjoin(missingVars, ', '));

Sum = Sum(:, needVars);

%% Load HBI (cbm)
assert(exist(hbi_file,'file')==2, 'Cannot find %s. Please run cbm\_hbi.', hbi_file);
S = load(hbi_file);
assert(isfield(S,'cbm'), 'No cbm variable in %s.', hbi_file);
cbm = S.cbm;

assert(isfield(cbm,'output') && isfield(cbm.output,'responsibility'), 'cbm.output.responsibility does not exist.');
resp = cbm.output.responsibility;

% Names of the models in the same order as cbm.output.responsibility
% Mentalizing, Level 0, Level 1 (diversity), Level 2 (knowledgeability) and Baseline
strategyNames = {'Mentalizing','Level0','DiversityHeuristic','ColourModel','Baseline'};
nModels_expected = numel(strategyNames);

assert(size(resp,2) == nModels_expected, ...
    'Responsibility dimension mismatch: expected N×%d but got %s.', ...
    nModels_expected, mat2str(size(resp)));

N_resp = size(resp,1);
N_sum  = height(Sum);

assert(N_resp == N_sum, ...
    'Subject count mismatch: Excel=%d rows, responsibility=%d rows.\n', ...
     N_sum, N_resp);

%% Cluster participants into strategies
[respMax, stratIdx] = max(resp, [], 2);

% Force column vectors
respMax  = respMax(:);
stratIdx = stratIdx(:);

% strategy categorical
strategy = categorical(strategyNames(stratIdx), strategyNames);
strategy = strategy(:);

% strategyU (with Uncertain)
strategyU = categorical(strategyNames(stratIdx), [strategyNames, {'Uncertain'}]);
strategyU(respMax < confThr) = categorical({'Uncertain'});
strategyU = strategyU(:);

StratTbl = table(strategy, strategyU, respMax, stratIdx, ...
    'VariableNames', {'strategy','strategyU','respMax','stratIdx'});

% Convert responsibility matrix into a table with meaningful column names
RespTbl = array2table(resp, 'VariableNames', ...
    {'r_mentalizing','r_level0','r_diversity','r_colour','r_baseline'});

A = [Sum, StratTbl, RespTbl];

% Remove rows with missing metrics
A = rmmissing(A, 'DataVariables', needVars);

% exclude Uncertain

A = A(A.strategyU ~= categorical({'Uncertain'}), :);


% Drop empty categories to avoid "jumping indices" in old functions
A.strategyU = removecats(A.strategyU);

%% Mixed Effect ANOVA
Y = A(:, {'div_t','nodiv_t'});

Between = table(categorical(A.strategy), 'VariableNames', {'strategyU'});

T = [Between Y];

T.strategyU = removecats(categorical(T.strategyU));


%% 

within = table(categorical(["div"; "nodiv"]), 'VariableNames', {'Condition'});

rm = fitrm(T, 'div_t,nodiv_t ~ strategyU', 'WithinDesign', within);

ranovatbl = ranova(rm, 'WithinModel', 'Condition');

betweentbl = anova(rm);

disp("=== Within-subject (Condition) + Interaction (Condition*strategyU) ===")
disp(ranovatbl)

disp("=== Between-subject main effect (strategyU) ===")
disp(betweentbl)


%% ANOVA
metrics = {'av_pdiv','av_nodiv','time','t','median_t','div_t','nodiv_t'};
minGroupN = 5;  % or your value

R = cell(0,8);                % {metric,N,K,F,pA,eta2,pLev,pKW}
PosthocAll = table();

%% 
for k = 1:numel(metrics)
    yName = metrics{k};
    y = A.(yName);
    g = A.strategyU;

    ok = ~isnan(y) & ~isundefined(g);
    y2 = y(ok);
    g2 = removecats(g(ok));

    [gID, gNames] = findgroups(g2);
    K = numel(gNames);
    if K < 2
        warning('%s: <2 groups; skipped.', yName); 
        continue;
    end

    grpN = accumarray(gID, 1);
    if any(grpN < minGroupN)
        warning('%s: group N < %d (%s).', yName, minGroupN, mat2str(grpN(:)'));
    end

    % Levene (median abs dev) via one-way ANOVA on deviations
    med  = splitapply(@median, y2, gID);
    pLev = anova1(abs(y2 - med(gID)), gID, 'off');

    % One-way ANOVA
    [pA, tbl, stats] = anova1(y2, gID, 'off');
    eta2 = tbl{2,2} / tbl{4,2};
    Fstat = tbl{2,5};

    % Kruskal-Wallis
    pKW = kruskalwallis(y2, gID, 'off');

    % Tukey-Kramer posthoc
    mc = multcompare(stats, 'CType','tukey-kramer', 'Display','off');
    gLabels = cellstr(categorical(gNames));   % robust conversion

    ph = table( ...
        repmat({yName}, size(mc,1), 1), ...
        gLabels(mc(:,1)), gLabels(mc(:,2)), ...
        mc(:,4), mc(:,3), mc(:,5), mc(:,6), ...
        'VariableNames', {'metric','group1','group2','diff','CI_low','CI_high','p'} );

    PosthocAll = [PosthocAll; ph];

    R(end+1,:) = {yName, numel(y2), K, Fstat, pA, eta2, pLev, pKW}; 
end

Results = cell2table(R, 'VariableNames', {'metric','N','kGroups','F', 'p_anova','eta2','p_levene','p_kw'});

disp('===ANOVA Summary===');
disp(Results);

disp('====Posthoc (Tukey-Kramer)===');
disp(PosthocAll);

%% descriptives result: mean/SD/SEM 

mean_nan = @(x) mean(x(~isnan(x)));
std_nan  = @(x) std(x(~isnan(x)));
n_nan    = @(x) sum(~isnan(x));

g = removecats(A.strategyU);
[Gid, gNames] = findgroups(g);
K = numel(gNames);

if iscategorical(gNames)
    gLabels = cellstr(gNames);
else
    gLabels = cellstr(categorical(gNames));
end

Desc = table(gLabels(:), accumarray(Gid,1), 'VariableNames', {'strategyU','N'});

for i = 1:numel(needVars)
    yName = needVars{i};
    y = A.(yName);

    m  = splitapply(mean_nan, y, Gid);
    sd = splitapply(std_nan,  y, Gid);
    n  = splitapply(n_nan,    y, Gid);
    sem = sd ./ sqrt(max(n,1));

    Desc.(['mean_' yName]) = m;
    Desc.(['sd_'   yName]) = sd;
    Desc.(['sem_'  yName]) = sem;
    Desc.(['n_'    yName]) = n;
end

disp('=== Group descriptives (mean, SEM)===');
disp(Desc);


