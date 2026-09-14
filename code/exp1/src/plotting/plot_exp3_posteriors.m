%% Exp3 posterior weights: eight square bar-chart panels
% User-set prediction parameters on the positive scale (AFTER exp).
% ce is retained for consistency with the models; it does not change w.
clear; clc;
paramValues.mentalizing = 5;         % [ce]
paramValues.vt = [5 1];              % [ce, n]
paramValues.dw = [5 1 5];           % [ce, n, e]
paramValues.kw = [5 1];             % [ce, n]
trialsToPlot = 1:8;                 % Alternatively, 9:16
saturationScale = 0.82;             % 1 = supplied colors; lower = less saturated
saveFigure = true;
projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
outputDir = fullfile(projectRoot, 'figures', 'exp3_posteriors');

names = {'mentalizing','vt','dw','kw'};
labels = {'M','VT','DW','KW'};
counts = [1 2 3 2];
for k = 1:4
    validateattributes(paramValues.(names{k}), {'numeric'}, ...
        {'row','real','finite','positive','numel',counts(k)});
end
validateattributes(trialsToPlot, {'numeric'}, ...
    {'row','numel',8,'integer','>=',1,'<=',16});
assert(numel(unique(trialsToPlot))==8, 'Select eight distinct trials.');
validateattributes(saturationScale, {'numeric'}, {'scalar','real','finite','>=',0,'<=',1});

%% Formulas from /Users/foxxis/Desktop/HBI/exp3/src/models
% Keep all 16 entries: KW swaps its conflict-trial weights in the second block.
posteriorWeight = zeros(16,4);
m = [0 0 1 1 .2735507 .7765487 .4667553 .5992063];
posteriorWeight(:,1) = [m m]';
n = paramValues.vt(2);
d = 2+2*n;
v = [1/d (1+2*n)/d .5 .5 1/d (1+2*n)/d 1/d (1+2*n)/d];
posteriorWeight(:,2) = [v v]';
n = paramValues.dw(2); e = paramValues.dw(3);
d = 2+2*n; de = 2+2*n*e;
v = [1/de (1+2*n*e)/de .5 .5 1/d (1+2*n)/d 1/d (1+2*n)/d];
posteriorWeight(:,3) = [v v]';
n = paramValues.kw(2);
eargain = 1.675676; footgain = 1.027027;
d = 2+(eargain+footgain)*n;
v = [1/d (1+(eargain+footgain)*n)/d ...
     (1+footgain*n)/d (1+eargain*n)/d ...
     1/(2+2*eargain*n) (1+2*eargain*n)/(2+2*eargain*n) ...
     1/(2+2*footgain*n) (1+2*footgain*n)/(2+2*footgain*n)];
posteriorWeight(:,4) = [v v([1 2 4 3 5 6 7 8])]';
assert(all(isfinite(posteriorWeight(:))) && ...
    all(posteriorWeight(:)>=0 & posteriorWeight(:)<=1));

%% Preserve the supplied hues and brightness while reducing HSV saturation
baseColors = [ ...
    .070 .465 .710; ...
    .830 .350 .140; ...
    .910 .570 .170; ...
    .850 .685 .140];
hsvColors = rgb2hsv(baseColors);
hsvColors(:,2) = saturationScale*hsvColors(:,2);
colors = hsv2rgb(hsvColors);

fig = figure('Color','w','Name','Experiment 3 posterior predictions', ...
    'Position',[100 40 680 1000]);
set(fig,'DefaultAxesFontName','Arial','DefaultTextFontName','Arial');
layout = tiledlayout(fig,4,2,'TileSpacing','compact','Padding','compact');
title(layout,'Experiment 3 | Model predictions','FontName','Arial', ...
    'FontSize',14,'FontWeight','bold');
% MATLAB tiles are row-major; use column-wise trial order as in the reference.
panelOrder = [1 5 2 6 3 7 4 8];
for panel = 1:8
    trial = trialsToPlot(panelOrder(panel));
    ax = nexttile(layout);
    bar(ax,1:4,posteriorWeight(trial,:),.72,'FaceColor','flat', ...
        'CData',colors,'EdgeColor',[.25 .25 .25],'LineWidth',.7);
    set(ax,'FontName','Arial','FontSize',9,'Box','off','TickDir','out', ...
        'LineWidth',.7,'XColor',[.2 .2 .2],'YColor',[.2 .2 .2], ...
        'XGrid','off','YGrid','on','GridColor',[.85 .85 .85],'GridAlpha',.35);
    xlim(ax,[.4 4.6]); ylim(ax,[0 1]);
    xticks(ax,1:4); xticklabels(ax,labels); yticks(ax,[0 .5 1]);
    ylabel(ax,'P(Type II | s_1, s_2)','Interpreter','tex', ...
        'FontName','Arial','FontSize',9);
    title(ax,sprintf('Trial %d',trial),'FontName','Arial','FontSize',10);
    pbaspect(ax,[1 1 1]);
end
if saveFigure
    if ~exist(outputDir,'dir'), mkdir(outputDir); end
    stem = sprintf('exp3_trials_%d_to_%d',trialsToPlot(1),trialsToPlot(end));
    exportgraphics(fig,fullfile(outputDir,[stem '.png']),'Resolution',300);
    exportgraphics(fig,fullfile(outputDir,[stem '.pdf']),'ContentType','vector');
    savefig(fig,fullfile(outputDir,[stem '.fig']));
    save(fullfile(outputDir,[stem '.mat']),'paramValues','posteriorWeight', ...
        'trialsToPlot','colors','saturationScale');
end
