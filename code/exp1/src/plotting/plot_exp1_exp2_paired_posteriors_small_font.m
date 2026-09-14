%% Exp1 and Exp2: posterior probabilities and participant contrasts (4 x 4)
% All parameter values below are AFTER exp. ce does not affect these weights.
clear; clc;
paramValues(1).mentalizing = 5;
paramValues(1).vt = [5 2];
paramValues(1).dw = [5 2 5];
paramValues(1).kw = [5 2];
% Exp2 defaults preserve the values in exp2/Posterior_plot.m.
paramValues(2).mentalizing = 5;
paramValues(2).vt = [5 2];
paramValues(2).dw = [5 2 5];
paramValues(2).kw = [5 2];
complementTrials = {[1 2], []};
% Display order in terms of original trial IDs; swap every Exp2 pair.
trialOrder = {1:8, [2 1 4 3 6 5 8 7]};
behaviorFiles = { ...
    '/Users/foxxis/Documents/social induction/experiment data/exp1/output100.xlsx', ...
    '/Users/foxxis/Documents/social induction/experiment data/exp2/final_data.csv'};
behaviorColumns = {{'d1','d2','d3','d4','d5','d6','d7','d8'}, ...
                   {'p1','p2','p3','p4','p5','p6','p7','p8'}};
saveFigure = true;
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
outputDir = fullfile(root, 'figures', 'exp1_exp2_paired');

%% Compute each experiment independently, using its own model definitions
results = struct([]);
allContrasts = [];
for experiment = 1:2
    names = {'mentalizing','vt','dw','kw'};
    counts = [1 2 3 2];
    for k = 1:4
        validateattributes(paramValues(experiment).(names{k}), {'numeric'}, ...
            {'real','row','numel',counts(k),'finite','positive'});
    end
    flipped = complementTrials{experiment};
    order = trialOrder{experiment};
    assert(isequal(sort(order),1:8), 'trialOrder must be a permutation of 1:8.');
    if ~isempty(flipped)
        validateattributes(flipped, {'numeric'}, {'vector','integer','>=',1,'<=',8});
    end
    rawWeights = modelWeights(experiment, paramValues(experiment));
    assert(all(isfinite(rawWeights(:))) && all(rawWeights(:)>=0 & rawWeights(:)<=1));
    weights = rawWeights;
    weights(flipped,:) = 1 - weights(flipped,:);
    weights = weights(order,:);
    contrasts = weights(1:2:8,:) - weights(2:2:8,:);
    assert(isfile(behaviorFiles{experiment}), 'Missing data: %s', behaviorFiles{experiment});
    tbl = readtable(behaviorFiles{experiment});
    columns = behaviorColumns{experiment};
    assert(all(ismember(columns,tbl.Properties.VariableNames)), 'Missing trial columns.');
    ratings = table2array(tbl(:,columns));
    if ~isnumeric(ratings), ratings = str2double(string(ratings)); end
    validateattributes(ratings, {'numeric'}, {'real','2d','ncols',8,'nonempty'});
    valid = isfinite(ratings) & ratings>=1 & ratings<=6 & ratings==fix(ratings);
    if any(~valid(:))
        warning('Exp%d: excluding %d invalid or missing ratings.',experiment,nnz(~valid));
    end
    probabilities = (ratings-0.5)/6;
    probabilities(~valid) = NaN;
    probabilities(:,flipped) = 1-probabilities(:,flipped);
    probabilities = probabilities(:,order);
    observed = probabilities(:,1:2:8)-probabilities(:,2:2:8);
    results(experiment).rawWeights = rawWeights;
    results(experiment).trialOrder = order;
    results(experiment).weights = weights;
    results(experiment).contrasts = contrasts;
    results(experiment).rawRatings = ratings;
    results(experiment).ratings = ratings(:,order);
    results(experiment).probabilities = probabilities;
    results(experiment).observed = observed;
    results(experiment).pairCounts = sum(isfinite(observed),1);
    fprintf('Exp%d valid counts: %d %d %d %d\n',experiment,results(experiment).pairCounts);
    displayed = contrasts(:,1:2);
    allContrasts = [allContrasts; displayed(:); observed(isfinite(observed))]; %#ok<AGROW>
end
deltaLimit = max(1,ceil((max(abs(allContrasts))+0.16)*10)/10);

%% Four columns: Exp1 posterior / contrast / Exp2 posterior / contrast
colors = [195 154 56;67 134 117;162 127 169;88 127 165]/255;
labels = {'M','VT','DW','KW'};
fig = figure('Color','w','Name','Experiments 1 and 2: paired posteriors', ...
    'Position',[40 40 1580 1050]);
set(fig,'DefaultAxesFontName','Arial','DefaultTextFontName','Arial');
heading(fig,[.055 .95 .90 .035],'Model predictions and participant behavior',16);
heading(fig,[.055 .89 .42 .03],'Experiment 1',12);
heading(fig,[.55 .89 .42 .03],'Experiment 2',12);
left = [.055 .30 .55 .795];
width = [.19 .175 .19 .175];
columnTitles = {'a   Posterior probability','b   Within-pair contrast', ...
                'c   Posterior probability','d   Within-pair contrast'};
for column = 1:4
    heading(fig,[left(column) .843 width(column)+.025 .03],columnTitles{column},10);
end
annotation(fig,'line',[.507 .507],[.105 .923],'Color',[.82 .82 .82],'LineWidth',.6);
for pair = 1:4
    trials = [2*pair-1 2*pair];
    bottom = .695-(pair-1)*.19;
    for experiment = 1:2
        column = 2*experiment-1;
        r = results(experiment);
        sourceTrials = r.trialOrder(trials);
        ax = axes(fig,'Position',[left(column) bottom width(column) .12]);
        hold(ax,'on');
        for k = 1:4
            bar(ax,k-.17,r.weights(trials(1),k),.30,'FaceColor',colors(k,:), ...
                'EdgeColor',colors(k,:),'LineWidth',1.1);
            bar(ax,k+.17,r.weights(trials(2),k),.30,'FaceColor','w', ...
                'EdgeColor',colors(k,:),'LineWidth',1.1);
        end
        styleAxes(ax,labels); ylim(ax,[0 1]); yticks(ax,[0 .5 1]);
        flipped = ismember(sourceTrials,complementTrials{experiment});
        if all(flipped)
            probabilityLabel = '1-P(TypeII|s_1,s_2)';
            deltaLabel = '\Delta(1-P)';
        elseif ~any(flipped)
            probabilityLabel = 'P(TypeII|s_1,s_2)';
            deltaLabel = '\Delta P';
        else
            probabilityLabel = 'Displayed probability'; deltaLabel = '\Delta probability';
        end
        ylabel(ax,probabilityLabel,'Interpreter','tex','FontName','Arial','FontSize',9);
        title(ax,sprintf('Pair %d | Trial %d vs %d',pair,sourceTrials),'FontSize',9);
        if pair==1 && experiment==1
            h1 = patch(ax,NaN,NaN,[.45 .45 .45],'EdgeColor',[.35 .35 .35]);
            h2 = patch(ax,NaN,NaN,'w','EdgeColor',[.35 .35 .35]);
            lg = legend(ax,[h1 h2], ...
                {'with-diversity condition','no-diversity condition'}, ...
                'Orientation','vertical','Box','off','FontSize',9, ...
                'FontName','Arial','Interpreter','none','AutoUpdate','off');
            set(lg,'Units','normalized','Position',[.75 .035 .22 .055]);
        end
        column = column+1;
        ax = axes(fig,'Position',[left(column) bottom width(column) .12]);
        hold(ax,'on'); yline(ax,0,'-','Color',[.35 .35 .35],'LineWidth',.8);
        for k = 1:2
            value = r.contrasts(pair,k);
            plot(ax,[k k],[0 value],'-','Color',colors(k,:),'LineWidth',2);
            plot(ax,k,value,'o','MarkerSize',6,'MarkerFaceColor',colors(k,:), ...
                'MarkerEdgeColor',colors(k,:));
            if abs(value)<.005, txt='0.00'; else, txt=sprintf('%+.2f',value); end
            if value>=-1e-6, offset=.07; align='bottom'; else, offset=-.07; align='top'; end
            text(ax,k,value+offset,txt,'HorizontalAlignment','center', ...
                'VerticalAlignment',align,'FontName','Arial','FontSize',8);
        end
        observed = r.observed(:,pair);
        observed = observed(isfinite(observed));
        x = 3*ones(size(observed));
        [~,~,groups] = unique(round(observed*6));
        for group = 1:max([0;groups])
            members = find(groups==group);
            if numel(members)>1
                spread = min(.24,.018*(numel(members)-1));
                x(members) = 3+linspace(-spread,spread,numel(members))';
            end
        end
        scatter(ax,x,observed,12,[.25 .25 .25],'filled', ...
            'MarkerFaceAlpha',.45,'MarkerEdgeColor','none');
        xline(ax,2.5,':','Color',[.75 .75 .75]);
        styleAxes(ax,{'M','VT','Data'});
        ylim(ax,[-deltaLimit deltaLimit]); yticks(ax,[-1 -.5 0 .5 1]);
        ylabel(ax,deltaLabel,'Interpreter','tex','FontName','Arial','FontSize',9);
        title(ax,sprintf('Trial %d - Trial %d',sourceTrials),'FontSize',9,'FontWeight','normal');
    end
end

if saveFigure
    if ~exist(outputDir,'dir'), mkdir(outputDir); end
    exportgraphics(fig,fullfile(outputDir,'exp1_exp2_paired.png'),'Resolution',300);
    exportgraphics(fig,fullfile(outputDir,'exp1_exp2_paired.pdf'),'ContentType','vector');
    savefig(fig,fullfile(outputDir,'exp1_exp2_paired.fig'));
    save(fullfile(outputDir,'exp1_exp2_paired.mat'),'results','paramValues', ...
        'complementTrials','trialOrder','behaviorFiles','behaviorColumns');
end

function w = modelWeights(experiment,p)
% Explicit experiment-specific weight formulas from src/models/*_model.m.
w = zeros(8,4);
if experiment==1
    w(:,1) = [.8;.2380952;.0625;.4571429;.7999999;.8;.7999999;.8];
else
    w(:,1) = [.6359223;.8769231;.3858696;.7142857;.5409836;.5526316;.2807018;.3636364];
end
n=p.vt(2); d=2+2*n; a=(1+2*n)/d;
if experiment==1, w(:,2)=[.5;1/d;1/d;1/d;a;a;.5;a];
else, w(:,2)=[a;a;1/d;.5;a;a;1/d;.5]; end
n=p.dw(2); e=p.dw(3); d=2+2*n; de=2+2*n*e; a=(1+2*n)/d; ae=(1+2*n*e)/de;
if experiment==1, w(:,3)=[.5;1/d;1/de;1/d;ae;a;.5;a];
else, w(:,3)=[a;ae;1/d;.5;a;ae;1/d;.5]; end
n=p.kw(2);
if experiment==1
    b=[1.870504 1.108434 1.108434 1.108434];
    y=[3.085106 1.870504 1.108434 1.870504];
    w(:,4)=[(1+y(1)*n)/(2+(b(1)+y(1))*n);1/(2+2*b(1)*n); ...
        1/(2+(b(2)+y(2))*n);1/(2+2*b(2)*n); ...
        (1+(b(3)+y(3))*n)/(2+(b(3)+y(3))*n);(1+2*b(3)*n)/(2+2*b(3)*n); ...
        (1+b(4)*n)/(2+(b(4)+y(4))*n);(1+2*b(4)*n)/(2+2*b(4)*n)];
else
    b=[1.604278 1.604278 1.410604 1.410604 1.040060 1.040060 1.463840 1.463840];
    y=b;
    w(:,4)=[(1+2*b(1)*n)/(2+2*b(1)*n);(1+(b(2)+y(2))*n)/(2+(b(2)+y(2))*n); ...
        1/(2+2*y(3)*n);(1+b(4)*n)/(2+(y(4)+b(4))*n); ...
        (1+2*b(5)*n)/(2+2*b(5)*n);(1+(b(6)+y(6))*n)/(2+(b(6)+y(6))*n); ...
        1/(2+2*b(7)*n);(1+y(8)*n)/(2+(y(8)+b(8))*n)];
end
end

function heading(fig,position,txt,size)
annotation(fig,'textbox',position,'String',txt,'EdgeColor','none', ...
    'FontName','Arial','FontSize',size,'FontWeight','bold');
end

function styleAxes(ax,labels)
set(ax,'FontName','Arial','FontSize',9,'Color','w','Box','off','LineWidth',.65, ...
    'TickDir','out','XColor',[.3 .3 .3],'YColor',[.3 .3 .3], ...
    'XGrid','off','YGrid','on','GridColor',[.85 .85 .85],'GridAlpha',.4);
xlim(ax,[.4 numel(labels)+.6]); xticks(ax,1:numel(labels)); xticklabels(ax,labels);
end
