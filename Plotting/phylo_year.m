clear 
addpath ~/GitHub/Plankton_IBM/Functions
bioinf_fcns  = bioinf_functions; 

cd ~/GitHub/Plankton_IBM/output
load('workspace.mat','eco_params','env_forcing')
phylodata = load('last_year_tree.mat');

srcnode = phylodata.srcnode2;
snknode = phylodata.snknode2;

%% create adjacency matrix and graph for sample tree
nadj = max(snknode); % 
adj  = sparse(srcnode,snknode,1,nadj,nadj);

% generate directed graph object of sampled tree
Graph_directed = digraph(adj);

% get index of terminal leaf nodes of live cells (i.e. nodes that are sinks but not sources) 
ilive = find(outdegree(Graph_directed)==0); % terminal node outdegree == 0 
%% sample terminal leaf nodes 
nsample = 250; % max nsample terminal leaf nodes
nsample = min(nsample,numel(ilive)); % decrease if not enough nodes
isample = sort(randsample(numel(ilive),nsample)); % sample without replacement

%% use srcnodes to backtrack descent of each extant agent to source node

% create cell array with cell for each sampled agent (terminal leaf node)
nodeid = cell(nsample,1);
for i=isample'
    nodeid{i}=snknode(ilive(i)); % index of terminal leaf node
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

%% extract data from structures
srcnode  = phylodata.srcnode2(nodes);
snknode  = phylodata.snknode2(nodes);
 divtime = phylodata.divtime(nodes);
   t_opt = phylodata.t_opt(nodes);
 rgbgene = phylodata.rgbgene(nodes,:);
deadtime = phylodata.deadtime(nodes);

% clear phylodata  

%% reindex src and snknode IDs to make snk nodes consectutive
% make snk node indices consecutive
[~,snknode2] = sort(snknode);
% reindex src nodes to match snk nodes
ii = sparse(size(snknode2,2),1); % create blank sparse matrix
ii(snknode,1) = snknode2;        % map original snk index to new snk index
srcnode2 = full(ii(srcnode));    % convert src nodes using snknode mapping

%% calculate generation time (birth time - birth time of parent)
gentime = divtime(snknode2)-divtime(srcnode2);
% normalise rgb gene (centered on origin)
maxclr  = max(abs(rgbgene));
minclr  = -maxclr;
plotclr = (rgbgene-minclr)./(maxclr-minclr);

%% create adjacency matrix and graph for sample tree
nadj = max(snknode2); % 
adj  = sparse(srcnode2,snknode2,gentime,nadj,nadj);

% generate directed graph object of sampled tree
Graph_directed = digraph(adj);
% create undirected graph (to calculate distances)
Graph_undirected = graph(adj+adj');

% get index of terminal leaf nodes of live cells (i.e. nodes that are sinks but not sources) 
ilive = find(outdegree(Graph_directed)==0); % terminal node outdegree == 0 

% calculate node divergence by time and by n generations
% time divergence (uses edge weights)
d    = distances(Graph_undirected,ilive,ilive);  
% generation divergence (ignores edge weights)
dgen = distances(Graph_undirected,ilive,ilive,'Method','unweighted'); 

%% draw color-coded trees of exact graph-based descent

radial = false; % draw trees using linear or polar projection

figure(101)
clf

subplot(211)
hold on
% plot sampled tree
h = plot(Graph_directed,'Layout','layered');
if radial
    rho=divtime(nodes)';
    theta=normalize(h.XData,'range').*1.975.*pi;
    [x,y] = pol2cart(theta,rho);
    axis square
    h.XData = x;
    axis off
else
    y=divtime;
    h.YData = y;
    axis ij
end
h.LineWidth=2;
h.NodeColor='none';
h.ArrowSize=0;
h.EdgeAlpha=1;
h.NodeLabel=[];
h.EdgeColor=plotclr(Graph_directed.Edges.EndNodes(:,2),:);
hold on
scatter(h.XData(ilive),h.YData(ilive),25,t_opt(ilive,:),'filled');
nyears=env_forcing.tmax./env_forcing.daysperyear;
set(gca,'XTick',[]);
set(gca,'YTick',linspace(0,env_forcing.tmax,nyears+1));
set(gca,'YTickLabel',linspace(0,nyears,nyears+1));
ylabel('Time (years)')
xlim([-2 nsample+2])
box on
cx=caxis;
set(gcf,'Color','w')


subplot(212)
% plot sampled tree
h = plot(Graph_directed,'Layout','layered');
if radial
    rho=divtime(nodes)';
    theta=normalize(h.XData,'range').*1.975.*pi;
    [x,y] = pol2cart(theta,rho);
    axis square
    h.XData = x;
    axis off
else
    y=divtime;
    h.YData = y;
    axis ij
end
h.LineWidth=2;
h.NodeColor='none';
h.ArrowSize=0;
h.EdgeAlpha=1;
h.NodeLabel=[];

crng = linspace(cx(1),cx(2),65);
t_edges=t_opt(Graph_directed.Edges.EndNodes(:,2));
t_edges(t_edges<cx(1))=cx(1);
t_edges(t_edges>cx(2))=cx(2);
cbin = discretize(t_edges,crng);
cmap=parula(64);
h.EdgeColor=cmap(cbin,:);
hold on
scatter(h.XData(ilive),h.YData(ilive),25,plotclr(ilive,:),'filled')
nyears=env_forcing.tmax./env_forcing.daysperyear;
set(gca,'XTick',[]);
set(gca,'YTick',linspace(0,env_forcing.tmax,nyears+1));
set(gca,'YTickLabel',linspace(0,nyears,nyears+1));
ylabel('Time (years)')
xlim([-2 nsample+2])
box on
set(gcf,'Color','w')
%%






