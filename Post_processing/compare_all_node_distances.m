clear 
addpath ~/GitHub/Plankton_IBM/Functions
bioinf_fcns  = bioinf_functions; 

cd ~/GitHub/Plankton_IBM/output/sinusoidal_10year_IBM_resting_verbose
load('workspace.mat','eco_params','env_forcing')

fntsz = 16;

% load('pruned_tree.mat');
load('last_year.mat');


srcnode  = phylogeny_table.srcnode;
snknode  = phylogeny_table.snknode;
%% sample terminal leaf nodes 
nsample = 50; % max nsample nodes
nsample = min(nsample,numel(snknode)); % decrease if not enough nodes
isample = sort(randsample(snknode,nsample)); % sample without replacement
% isample gives location of extant agents in phylogeny_index2 and matching rgb and genome darrays


%% use srcnodes to backtrack descent of each extant agent to source node

% create cell array with cell for each sampled agent (terminal leaf node)
nodeid = cell(nsample,1);
for i=1:nsample
    nodeid{i}=isample(i); % index of terminal leaf node
    % iterate back through generations
    while nodeid{i}(end)~=1 % until root is reached
        nodeid{i}(end+1)=srcnode(nodeid{i}(end));
    end
end
% find unique nodes among all paths back from sampled terminal nodes
[nodes,~] = unique(horzcat(nodeid{:}));
% this is indexes all the nodes in the phylogenetic tree

clear nodeid 
%%

phylogeny_table = phylogeny_table(nodes,:);

%% reindex src and snknode IDs to make snk nodes consectutive
srcnode = phylogeny_table.srcnode;
snknode = phylogeny_table.snknode;
% make snk node indices consecutive
[~,snknode2] = sort(snknode);
% reindex src nodes to match snk nodes
ii = sparse(size(snknode2,2),1); % create blank sparse matrix
ii(snknode,1) = snknode2;        % map original snk index to new snk index
srcnode2 = full(ii(srcnode));    % convert src nodes using snknode mapping
phylogeny_table.srcnode = srcnode2;
phylogeny_table.snknode = snknode2;

%% extract variables from table

srcnode  = phylogeny_table.srcnode;
snknode  = phylogeny_table.snknode;
 divtime = phylogeny_table.divtime;
gentime  = divtime(snknode)-divtime(srcnode);
deadtime = phylogeny_table.deadtime;
   t_opt = phylogeny_table.t_opt;
 rgbgene = phylogeny_table.rgb_genome;

bingene_exists = ismember('binary_genome', phylogeny_table.Properties.VariableNames);
if bingene_exists
     bingene = phylogeny_table.binary_genome;
end

%%  normalise rgb gene (centered on origin)
maxclr  = max(abs(rgbgene));
minclr  = -maxclr;
plotclr = (rgbgene-minclr)./(maxclr-minclr);

%% create adjacency matrix and graph for sample tree
nadj = max(snknode); % 
adj  = sparse(srcnode,snknode,gentime,nadj,nadj);

% generate directed graph object of sampled tree
Graph_directed = digraph(adj);
% create undirected graph (to calculate distances)
Graph_undirected = graph(adj+adj');

% calculate node divergence by time and by n generations
% time divergence (uses edge weights)
d    = distances(Graph_undirected);  
% generation divergence (ignores edge weights)
dgen = distances(Graph_undirected,'Method','unweighted'); 

if bingene_exists
    % estimate divergence from molecular clock
    % find live cell index matching sample index
    % convert condensed-format genome to full binary strings
    genome_string = bioinf_fcns.print_genomes(bingene);
    
    L = eco_params.ngenome.*64;
    
    % generate hamming-distance matrix
    p = pdist(genome_string-'0','Hamming');
    p = squareform(p); % make square matrix
    X = p.*L./ eco_params.pneutral; % convert from fraction of bits different to total bits different
                                    % also account for mutation probability
    % 2-base jukes-cantor correction
    pcorr=p;
    pcorr(pcorr>0.5)=NaN;
    jk = -1/2 * log(1 - pcorr * 2/1);
    jk = jk.*L; 
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig1=figure(1);
fig1.Position = [152 700 900 700];



if bingene_exists
    
    % compare exact distances with estimates
    sh=subplot(223);
    cla
    hold off
    scatter(dgen(:), X(:),15,'k','filled','MarkerFaceAlpha',0.1,'MarkerEdgeColor','none');
    hold on
    scatter(dgen(:),jk(:),15,'r','filled','MarkerFaceAlpha',0.1,'MarkerEdgeColor','none');

    p = linspace(0,max(dgen(:))./L,100);
    djk = -1/2 * log(1 - p * 2);
    E = sqrt(p.*(1-p)./(L*((1-2.*p).^2)));
    pl1=plot(djk.*L,p.*L,'-','LineWidth',2,'Color',[0.5 0.5 0.5]);
    plot((djk+E).*L,p.*L,':','LineWidth',2,'Color',[0.5 0.5 0.5]);
    plot((djk-E).*L,p.*L,':','LineWidth',2,'Color',[0.5 0.5 0.5]);

    pl2=plot(djk.*L,djk.*L,'-','LineWidth',2,'Color',[0.75 0 0]);
    plot(djk.*L,(djk+E).*L,':','LineWidth',2,'Color',[0.75 0 0]);
    plot(djk.*L,(djk-E).*L,':','LineWidth',2,'Color',[0.75 0 0]);
    
    xlim([0 max(dgen(:))]);
    xlabel('True divergence (generations)')
    ylabel('Estimated divergence (generations)')
    axis square
    box on
    title('(c) Estimated vs. true divergence')
    legend([pl2 pl1],'Corrected','Uncorrected','Location','SouthEast')

    
    subplot(224)
    cla
    hold off
    scatter(d(:)./365,X(:),15,'k','filled','MarkerFaceAlpha',0.1,'MarkerEdgeColor','none');
    hold on
    scatter(d(:)./365,jk(:),15,'r','filled','MarkerFaceAlpha',0.1,'MarkerEdgeColor','none');
    xlabel('True divergence  (years)')
    ylabel('Estimated divergence (generations)')
    axis square
    box on
    title('(c) Estimated vs. true (time) divergence')
    
    drawnow
end
