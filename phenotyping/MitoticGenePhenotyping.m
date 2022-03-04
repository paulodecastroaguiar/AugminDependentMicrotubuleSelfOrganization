% Genes with implications on mitosis – phenotypic similarities
% From
% "Augmin-dependent microtubule self-organization drives kinetochore fiber
% maturation in mammals"
% by
% Ana C. Almeida, Joana Oliveira, Danica Drpic, Liam P. Cheeseman, Joana
% Damas, Harris A. Lewin, Denis M. Larkin, Paulo Aguiar, António J.
% Pereira, Helder Maiato
% 
% maiato@i3s.up.pt


%% Clean up
clear; close all; clc


%% Choose fingerprint type
% parameters in fingerprint gets score 1 if the fraction of cells with
% reported events is > fingerprint_threshold (otherwise gets score 0)
% fingerprint_threshold = 0.25; 



%% Load data
filename = 'IM_live_cell_analysis_v3_AH.xlsx';

%conditions = sheetnames( filename ); 
[~,conditions,~] = xlsfinfo(filename);
conditions = convertCharsToStrings( conditions' );
conditions(1) = []; % discard "Glossary"

M = readtable( filename, 'Sheet', 'CONTROL' );
parameters = M.Properties.VariableNames;
parameters = parameters(5:end); % discard first 4 columns (CELLID, NEBD-MET, MET-AO, NEBD-AO)

% fill data structure with raw data
for c = 1:numel(conditions)
    M = readtable( filename, 'Sheet', conditions(c) );
    data(c).condition = conditions(c);
    data(c).all       = table2array( M(:,5:end) );
end


%% Show distributions
pooled_data = [];
for p = 1:numel(parameters)    
    tmp = [];
    for c = 1:numel(conditions)
        tmp = [tmp; data(c).all(:,p)];
    end
    pooled_data(:,p) = tmp;
end

figure
sgtitle('Global scores distributions for all cells');
for p = 1:numel(parameters)
    subplot(1,numel(parameters),p)
    histogram( pooled_data(:,p) )
    xlabel( '#events' )
    xticks([0,1])
    ylim([0,1400])
    title(parameters(p), 'Interpreter', 'none')
end


%% Produce fingerprints
for c = 1:numel(conditions)
    
    n_cells = size( data(c).all, 1 );
    data(c).fingerprint = zeros( 1, numel(parameters) );
    for p = 1:numel(parameters)
        values = data(c).all(:,p);
        ind = find( values > 0 );
        % fraction of cells with non-zero events
        fraction = numel(ind) / n_cells;
        data(c).fingerprint(p) = fraction;
%         if fraction >= fingerprint_threshold
%             data(c).fingerprint(p) = 1;
%         end            
    end

end


%% Fingerprints - bars
% CONTROL
figure
sgtitle('Phenotype fingerprint for Control');
c = 1;
hold on
bar( 1:numel(parameters), data(c).fingerprint );
xticks(1:numel(parameters))
xticklabels( parameters )
set(gca,'TickLabelInterpreter','none');
xtickangle( 45 );
ax = gca;
ax.XAxis.FontSize = 8;
xlim([0.5, 8.5])
hold off
ylabel( 'normalized' )
%title( conditions(c), 'Interpreter', 'none');

% ALL THE OTHER CONDITIONS
figure
sgtitle('Phenotype fingerprints for all conditions');
for c = 2:numel(conditions)
    subplot(6,11,c-1);    
    hold on
    bar( 1:numel(parameters), data(c).fingerprint );
    xticks(1:numel(parameters))
    xticklabels( parameters )
    set(gca,'TickLabelInterpreter','none');
    xtickangle( 45 );
    ax = gca;
    ax.XAxis.FontSize = 5;
    xlim([0.5, 8.5])
    hold off
    ylabel( 'normalized' )
    title( conditions(c), 'Interpreter', 'none');
end


%% Fingerprints - radar plots
edges = 0:2*pi/numel(parameters):2*pi;

% CONTROL
figure
c = 1;
polarhistogram( 'BinEdges', edges, 'BinCounts', data(c).fingerprint) 
thetaticks( 360/(2*pi)*( edges(1:end-1) + edges(2)/2 ) );
thetaticklabels( {'A','B','C','D','E','F','G','H'} );
rlim([0,1])
rticks([0.0, 0.5, 1.0])
title( conditions(c) );

% ALL THE OTHER CONDITIONS
figure
sgtitle('Phenotypes fingerprints for all conditions');
for c = 2:numel(conditions)
    subplot(6,11,c-1);
    polarhistogram( 'BinEdges', edges, 'BinCounts', data(c).fingerprint) 
    thetaticks( 360/(2*pi)*( edges(1:end-1) + edges(2)/2 ) );
    thetaticklabels( {'A','B','C','D','E','F','G','H'} );
    ax = gca;
    ax.ThetaAxis.FontSize = 5;
    txt = char( conditions(c) );
    rlim([0,1])
    rticks([0.0, 0.5, 1.0])
    title( txt(3:end), 'Interpreter','none');
end



%% PART 3 - PERFORM CLUSTER ANALYSIS
fingerprints_all = [];
for c = 1:numel(conditions)
    fingerprints_all = [ fingerprints_all; data(c).fingerprint ];
end

tree = linkage( fingerprints_all, 'average', 'euclidean');

labels = char( conditions );
labels = labels(:,3:end);
labels(1,1:7) = 'CONTROL';

figure
h = dendrogram(tree, 0, 'Labels', labels, 'Orientation', 'top', 'ColorThreshold', 0.45);
xtickangle( 45 );
ax = gca;
ax.XAxis.FontSize = 8;
set(gca,'TickLabelInterpreter','none');
set(h,'LineWidth',2)

