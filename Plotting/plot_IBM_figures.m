clear 
% close all

addpath Functions
bioinf_fcns  = bioinf_functions; 

% outdir = 'output';
% cd(outdir)
load('workspace.mat','eco_params','env_forcing')

fntsz = 14;

savefigs=true;

%    load('pruned_tree.mat');
    load('last_year.mat');

%%
srcnode = phylogeny_table.srcnode;
snknode = phylogeny_table.snknode;

phylogeny_index = find(phylogeny_table.TLN);

%% sample terminal leaf nodes 
nsample =100; % max nsample terminal leaf nodes
nsample = min(nsample,numel(phylogeny_index)); % decrease if not enough nodes
isample = sort(randsample(numel(phylogeny_index),nsample)); % sample without replacement
% isample gives location of extant agents in phylogeny_index2 and matching rgb and genome darrays

% extract sample of extant agents and corresponding rgb and genome data
phylogeny_terminal = phylogeny_index(isample);

%% use srcnodes to backtrack descent of each extant agent to source node

% create cell array with cell for each sampled agent (terminal leaf node)
nodeid = cell(nsample,1);
for i=1:nsample
    nodeid{i}=phylogeny_terminal(i); % index of terminal leaf node
    % iterate back through generations
    while nodeid{i}(end)~=1 % until root is reached
        nodeid{i}(end+1)=srcnode(nodeid{i}(end));
    end
end
% find unique nodes among all paths back from sampled terminal nodes
[nodes,~] = unique(horzcat(nodeid{:}));
% this is indexes all the nodes in the phylogenetic tree

% clear nodeid 
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

% get index of terminal leaf nodes of live cells
ilive = full(ii(phylogeny_terminal));
for i=1:size(nodeid,1)
    nodeid{i} = full(ii(nodeid{i}));
end
    


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
 

%% trait variance along lineages

IBM_output=load('IBM_output.mat');

T=linspace(0,max(IBM_output.tout),10000);
f999=figure(999);
clf
for i=1:size(nodeid,1)
    nds = nodeid{i}';
    temp = t_opt(nds);
    time = divtime(nds);
%     [f,xi] = ksdensity(temp);
%     plot(xi,f,'LineWidth',1,'Color',[0 0 0 0.01]);
    plot(time./360,temp,'LineWidth',1,'Color',[0 0 0 0.01])
    hold on
    
    p_interp(i,:) = interp1(time,temp,T,'linear',temp(1));
    
    
    
    lin_mean(i) = mean(t_opt(nodeid{i}));
    lin_stdv(i) = std(t_opt(nodeid{i}));
end

percentile=5;

plot(T./360,prctile(p_interp,percentile),'k-','LineWidth',1)
plot(T./360,prctile(p_interp,100-percentile),'k-','LineWidth',1)

% scatter(T,env_forcing.Tfunc(T),1,'r','filled','MarkerEdgeAlpha',0.1)

bio_out=IBM_output.yout(:,2:end-1);
cumperc=cumsum(100.*bio_out./sum(bio_out,2),2);
ycoord = linspace(eco_params.T_bins(1),eco_params.T_bins(end),size(cumperc,2));
contour(IBM_output.tout./360,ycoord,cumperc',[percentile 100-percentile],'LineColor',[1 1 1].*0.2,'LineWidth',1);
set(gca,'XTick',0:ceil(env_forcing.tmax./360))
ylim([12 18])
xlabel('Time (years)')
ylabel('Thermal optimum ^\circC')

set(gcf,'Color','w')
if savefigs
    exportgraphics(f999,'lineage_distribution.png','Resolution',450)
end


disp(['Mean of lineage means               = ' num2str(mean(lin_mean))])
disp(['Mean of lineage standard deviations = ' num2str(mean(lin_stdv))])
disp(['Standard deviation of lineage means = ' num2str(std(lin_mean))])

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
d    = distances(Graph_undirected,ilive,ilive)./2;  
% generation divergence (ignores edge weights)
dgen = distances(Graph_undirected,ilive,ilive,'Method','unweighted')./2; 

if bingene_exists
    % estimate divergence from molecular clock
    % find live cell index matching sample index
    % convert condensed-format genome to full binary strings
    genome_string = bioinf_fcns.print_genomes(bingene(ilive,:));
    
    L = eco_params.ngenome.*53;
    
    % generate hamming-distance matrix
    p = pdist(genome_string-'0','Hamming');
    p = squareform(p); % make square matrix
    X = p.*L./ eco_params.pneutral./2; % convert from fraction of bits different to total bits different
                                       % also account for mutation probability
    % 2-base jukes-cantor correction
    pcorr=p;
    pcorr(pcorr>0.5)=NaN;
    jk = -1/2 * log(1 - pcorr * 2/1);
    jk = jk./eco_params.pneutral;
end

% also calculate distance matrix from rgb gene
Xrgb = pdist(rgbgene(ilive,:),'euclidean');
Xrgb = squareform(Xrgb);


figure(333)
hist(p(:),100)
%

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%% draw color-coded trees of exact graph-based descent

radial = false; % draw trees using linear or polar projection

f101=figure(101);
f101.Position=[1428 584 730 761];
clf

% subplot(121)
hold on
% plot sampled tree
h1 = plot(Graph_directed,'Layout','layered','AssignLayers','alap');
if radial
    rho=divtime';
    theta=normalize(h1.XData,'range').*1.975.*pi;
    [x,y] = pol2cart(theta,rho);
    h1.XData = x;
    h1.YData = y;
    scatter(0,0,100,'k','filled')
    text(0,0,'  Initial seed','FontSize',fntsz+5)
    hold on
    scatter(h1.XData,h1.YData,50,'k');
    scatter(h1.XData,h1.YData,50,t_opt,'filled');
    f101.Position=[1428 998 900 900];
    axis square
%     caxis([14.75 15.25])
else
    h1.YData = h1.XData;
    h1.XData = divtime;
    [~,sortnodes] = sort(h1.YData(ilive));
    hold on
    axis ij
    scatter(h1.XData(ilive),h1.YData(ilive),25,'k');
    scatter(h1.XData(ilive),h1.YData(ilive),25,t_opt(ilive,:),'filled');
end
h1.LineWidth=2;
h1.NodeColor='none';
h1.ArrowSize=0;
h1.EdgeAlpha=1;
h1.NodeLabel=[];
h1.EdgeColor=plotclr(Graph_directed.Edges.EndNodes(:,2),:);
colormap(redblue)
if radial
    axis off
else
    nyears=env_forcing.tmax./env_forcing.daysperyear;
    set(gca,'YTick',[]);
    set(gca,'XTick',linspace(0,env_forcing.tmax,min(11,nyears+1)));
    set(gca,'XTickLabel',linspace(0,nyears,min(11,nyears+1)));
    xlabel('Time (years)')
    ylim([-10 nsample+10])
    title(['Known phylogeny for ' num2str(env_forcing.tmax/env_forcing.daysperyear) ' year run of IBM'])
    xlim([0 env_forcing.tmax])
%     caxis([0 30])
    cx=caxis;
    set(gcf,'Color','w')
    ch=colorbar('Location','eastoutside','FontSize',fntsz);
    set(get(ch,'ylabel'),'string','Temperature (^\circC)','FontSize',fntsz);
end
caxis([min(t_opt) max(t_opt)])
box on
pos=get(gca,'Position');
set(gca,'Position',pos);
set(gca,'FontSize',fntsz)



if savefigs
    if radial
        exportgraphics(f101,'radial.png','Resolution',450)
    else
        exportgraphics(f101,'phylogeny.png','Resolution',450)
    end
end

%%
IBM_output=load('IBM_output.mat');

if strcmp(func2str(env_forcing.Tfunc),'constant')
    ylimits = [13 17];
else
    ylimits = [8 22];
end 


bio_out=IBM_output.yout(:,2:end-1);
perc=(100.*bio_out./sum(bio_out,2));
ycoord = linspace(eco_params.T_bins(1),eco_params.T_bins(end),size(perc,2));


f102=figure(102);
f102.Position=[91 900 600 400];

clf
% imagesc(IBM_output.tout,ycoord,perc');
cmap=flipud(gray(8));
contourf(IBM_output.tout,ycoord,perc',[0:8],'LineColor',[1 1 1].*0.2,'LineWidth',0.25);
hold on
caxis([0 8])
colormap(cmap)

% draw descent tree
hold on
h2 = plot(Graph_directed,'Layout','layered','AssignLayers','alap');
h2.XData = divtime;
h2.YData = t_opt;
h2.LineWidth=1;
h2.NodeColor='none';
h2.ArrowSize=0;
h2.EdgeAlpha=1;
h2.NodeLabel=[];
h2.EdgeColor=plotclr(Graph_directed.Edges.EndNodes(:,2),:);




nyears=env_forcing.tmax./env_forcing.daysperyear;
axis([0 env_forcing.tmax ylimits])
set(gca,'YTick',ylimits(1):ylimits(end))
set(gca,'XTick',linspace(0,env_forcing.tmax,nyears+1));
set(gca,'XTickLabel',linspace(0,nyears,nyears+1));
xlabel('Time (years)')
box on
cx=caxis;
set(gcf,'Color','w')
ylabel('Temperature (^\circC)')
axis xy
set(gca,'FontSize',fntsz)
ch=colorbar('Location','west');
ch.Position(4) = ch.Position(4).*0.25;
ch.Ticks = 0:2:8;
ch.TickLabels{end} = [ch.TickLabels{end} '%'];


T=linspace(0,max(IBM_output.tout),10000);
scatter(T,env_forcing.Tfunc(T),1,'r','filled','MarkerEdgeAlpha',0.1)

drawnow
set(gcf,'Color','w')
if savefigs
    exportgraphics(f102,'time_series.png','Resolution',450)
end

%%
if ~strcmp(func2str(env_forcing.Tfunc),'squarewave')
    if strmatch(func2str(env_forcing.Tfunc),'constant')
        N = eco_params.gamma_0_mort.*eco_params.kN_half_sat./(eco_params.mu0-eco_params.gamma_0_mort);
%         N = min(IBM_output.yout(end-359:end,1))./eco_params.V;
    else
        N = IBM_output.yout(:,1)'./eco_params.V;
    end
    
    T   = env_forcing.Tfunc(1:env_forcing.tmax);
    T_phen = ycoord';
    
    F_T = exp( -((T - T_phen)./eco_params.theta_norm).^2 );
    F_N = Nutrient_function(N,eco_params);
    
    if false
        % Transport matrix
        eco_params.deltaM = eco_params.sigmaM.^2 ...      % genotype transfer function
            ./ ( 3.*eco_params.deltaT_opt.^2 );
        n      = eco_params.jpmax;                      % n 'genotypes'
        delta  = ones(n,1).*eco_params.delta;          % expanded vector of genotype transfer function
        TT     = spdiags([delta delta],[-1 1], n, n); % set off-diagonals to delta M
        TT     = spdiags(-sum(TT,2),0, TT);             % set diagonal to conserve mass
        TT = TT + speye(size(TT));
    end
    
    
    fitness =    (eco_params.mu0 .* F_T .* F_N) - eco_params.gamma_0_mort;
    fitness(fitness>0)=-eps;
    % fitness = TT*(eco_params.mu0 .* F_T .* F_N) - eco_params.gamma_0_mort;
    
    f654 = figure(654);
    f654.Position=[91 900 600 400];
    
    clf
    if strmatch(func2str(env_forcing.Tfunc),'constant')
        excl_time = log10([30 [1 10 100 1000].*360]);
        contourf(IBM_output.tout,ycoord,log10(-1./fitness),excl_time,'LineColor','none');
        caxis([min(excl_time) max(excl_time)])
        cmap=flipud(gray(5));
        colormap(cmap(1:end-1,:))
        ch=colorbar('Location','west');
        ch.Position(4) = ch.Position(4).*0.25;
        ch.Ticks = linspace(min(caxis),max(caxis),5);
        ch.TickLabels = {' ','Year','Decade','Century'};
        caxis([min(ch.Ticks) max(ch.Ticks)])
        hold on
        title('Ecological exclusion timescale')
    else
        imagesc(IBM_output.tout,ycoord,fitness);
        hold on
        contour(IBM_output.tout,ycoord,perc',[0:8],'LineColor',[1 1 1].*0.2,'LineWidth',0.25);
        colormap(redblue(64))
        ch=colorbar('Location','west');
        ch.Position(4) = ch.Position(4).*0.25;
        caxis([-0.1 0.1])
        title('Invasion fitness (1/day)')
    end
    
    set(gca,'XTick',linspace(0,env_forcing.tmax,nyears+1));
    set(gca,'XTickLabel',linspace(0,nyears,nyears+1));
    ylim(ylimits)
    
    
    
    hold on
    h2 = plot(Graph_directed,'Layout','layered','AssignLayers','alap');
    h2.XData = divtime;
    h2.YData = t_opt;
    h2.LineWidth=1;
    h2.NodeColor='none';
    h2.ArrowSize=0;
    h2.EdgeAlpha=1;
    h2.NodeLabel=[];
    h2.EdgeColor=plotclr(Graph_directed.Edges.EndNodes(:,2),:);
    T=linspace(0,max(IBM_output.tout),10000);
    % scatter(T,env_forcing.Tfunc(T),1,'k','filled','MarkerEdgeAlpha',0.1)
    axis xy
    
    xlim([0 env_forcing.tmax])
    if strmatch(func2str(env_forcing.Tfunc),'constant')
        contour(IBM_output.tout,ycoord,log10(-1./fitness),[1 1].*log10(360.*15),'k--','LineW',1);
    end
    
    
    ylabel('Temperature (^\circC)')
    set(gca,'YTick',ylimits(1):ylimits(end))
    
    set(gcf,'Color','w')
    if savefigs
        exportgraphics(f654,'fitness.png','Resolution',450)
    end
end
%% bioinformatics toolbox

clear bingenes

genmat=genome_string(sortnodes,:)-'0';
for i=1:nsample
    bingenes(i).Sequence=char(97+genmat(i,:)-'0');
    bingenes(i).Header=num2str(i);
end


UPGMAtree = seqlinkage(jk.*L,'UPGMA',bingenes);

leafnames   = strvcat(get(UPGMAtree,'LeafNames'));
sortnodes = str2num(leafnames);
[~,isort] = sort(sortnodes);
UPGMAtree = reorder(UPGMAtree,isort,'approximate',true);

h1 = plot(UPGMAtree,'orient','left');


% axis xy
title('UPGMA Distance Tree of agents using Jukes-Cantor model');
ylabel('Evolutionary distance')
if savefigs
    print -depsc -tiff -r300 -painters UPGMA.eps
end

% T = cluster(UPGMAtree,'MaxClust',5);

Xdist=[h1.LeafDots.XData h1.BranchDots.XData];
Ydist=[h1.LeafDots.YData h1.BranchDots.YData];
% Xdist=[h2.LeafDots.XData h2.BranchDots.XData];
% Ydist=[h2.LeafDots.YData h2.BranchDots.YData];
%%
f303=figure(303);
f303.Position = [139 103 1239 468];
clf

tree = UPGMAtree;


sh2=subplot(161);
adj = getmatrix(tree);
g2 = digraph(adj);
% add extra nodes to square up dendrogram
XDta = Xdist;
YDta = Ydist;
c = size(g2.Nodes,1); 
for i=1:size(g2.Nodes,1);
    preID = predecessors(g2,i);
    if numel(preID)>0
        c =c+1;
        g2 = addedge(g2,preID,c);
        XDta=[XDta XDta(preID)];
        YDta=[YDta YDta(i)];
        
        g2 = addedge(g2,c,i);
        g2 = rmedge(g2,preID,i);
    end
end
g=plot(g2);

g.XData=XDta;
g.YData=YDta;
% xlim([0 450])

g.LineWidth=2;
g.NodeColor='none';
g.ArrowSize=0;
g.EdgeAlpha=1;
g.NodeLabel=[];
g.EdgeColor='k';
axis ij
% [~,ind]=sort(outperm);
hold on
axis off
ylim([0.5 nsample+0.5])
title('(a) estimated phylogeny')

leafnames   = strvcat(get(tree,'LeafNames'));
sortnodes = flipud(str2num(leafnames));
cbin = discretize(t_opt(ilive(sortnodes)),linspace(min(t_opt),max(t_opt),64+1));

genmat=genome_string(sortnodes,:)-'0';
rgbmat=zeros(nsample,size(genmat,2),3);
for i=1:nsample
    clr=plotclr(ilive(sortnodes(i)),:)';
    rgbmat(i,:,:) = (clr*genmat(i,:))';
end
% rgbmat(rgbmat==0)=1;
subplot(1,6,[2 6])
imagesc(rgbmat)
axis xy
% imagesc(sign(rgbmat))
axis tight
box on
set(gca,'XTick',[],'YTick',[])
xlim([1 128])

cmap=redblue(64);
ch=colorbar('Location','eastoutside','Ticks',[],'FontSize',fntsz);
ylabel(ch, 'Thermal optimum')
colormap(cmap(cbin,:))
title('(b) neutral genome')

sh2.Position=[0.15  0.1100    0.1023    0.8150];
set(gcf,'Color','w')

if savefigs
    exportgraphics(f303,'genome.png','Resolution',450)
end



%%

nsample = 20; % max nsample terminal leaf nodes
iextant = find(IBM_output.y(2:end-1));
nextant = numel(iextant);
nsample = min(nsample,nextant); % decrease if not enough nodes
isample = sort(randsample(iextant,nsample)); % sample without replacement

f111=figure(111)
hold on
xx=sort(Xdist);
xx=xx-max(xx);
xx=xx(1:nsample);

plot(xx,1:nsample,'LineWidth',2)

thalf = 2.*nextant/nsample;
plot(xx,nsample./(1-xx./thalf),'k','LineWidth',2);

xlabel('Generations before end of simulation')
ylabel('Coalescences')

if savefigs
    export_fig -r450 genome.png
    print -depsc -tiff -r300 -painters .eps
    exportgraphics(f111,'coalesence.png','Resolution',450)
end

%%




dmat{1}=dgen(sortnodes,sortnodes);
dmat{2}=jk(sortnodes,sortnodes).*L./2;
titl={'(a) Known divergence','(b) Estimated divergence'};

fig1=figure(202);
fig1.Position = [152 700 900 700];
clf
cmap=redblue(64);

for i=1:2
    subplot(2,2,i)
    imagesc(dmat{i});
    axis square
    axis xy
    ch=colorbar('FontSize',fntsz);
    ylabel(ch, 'Divergence (generations)')
    ax1_pos = get(gca,'Position'); % position of first axes
    set(gca,'Position',ax1_pos,'XTick',[],'YTick',[]);
    caxis([0 max([dmat{1}(:);dmat{2}(:)])])
    cx=caxis;
    title(titl{i},'FontSize',fntsz)
    
    % rgb gene
    ax2 = axes('Position',ax1_pos,'XTick',[],'YTick',[]);
    axis square
    axis off
    ch=colorbar('Location','southoutside','Ticks',[],'FontSize',fntsz);
    ylabel(ch, 'RGB gene')
    colormap(ax2,plotclr(ilive(sortnodes),:))
    ax2.Position = ax1_pos;
    
    % thermal optima
    
    ax3 = axes('Position',ax1_pos,'XTick',[],'YTick',[]);
    axis square
    axis off
    ch=colorbar('Location','westoutside','Ticks',[],'FontSize',fntsz);
    ylabel(ch, 'Thermal optimum')
    colormap(ax3,cmap(cbin,:))
    ax3.Position = ax1_pos;
end

if bingene_exists
    
    % compare exact distances with estimates
    sh1=subplot(223);
    sh1.Position(1) = sh1.Position(1) - 0.03;
    scatter(dgen(:),X(:),15,'k','filled','MarkerFaceAlpha',0.1,'MarkerEdgeColor','none');
    hold on
    scatter(dgen(:),jk(:).*L./2,15,'r','filled','MarkerFaceAlpha',0.1,'MarkerEdgeColor','none');
    
    xi = linspace(0,max(dgen(:)),100); % true number of generations
    nonsat_ham=xi.*2./L;              % non-saturated Hamming distance (double distance)
    psat = 1/2.*(1 - exp(-2.*(nonsat_ham))); % estimated saturated  hamming distance
    
    ijk = linspace(0,i,100);
    xjk = psat;
    yjk = -1/2 * log(1 - xjk * 2/1);
    ejk = sqrt(xjk.*(1-xjk)./(L*((1-2.*xjk).^2)));
    
    pl1 =plot( xi,psat.*L./2,'-','LineWidth',2,'Color',[0.5 0.5 0.5]);

    pl2 =plot(xi, yjk     .*L./2,'-','LineWidth',2,'Color',[0.75 0 0]);
    pl2a=plot(xi,(yjk+ejk).*L./2,':','LineWidth',2,'Color',[0.75 0 0]);
    pl2b=plot(xi,(yjk-ejk).*L./2,':','LineWidth',2,'Color',[0.75 0 0]);
   
    xlim([0 max(dgen(:))]);
    xlabel('Known divergence (generations)')
    ylabel('Estimated divergence (generations)')
    axis square
    box on
    title('(c) Clock vs. generations')
    legend([pl2 pl1],'Corrected','Uncorrected','Location','SouthEast')
    set(gca,'FontSize',fntsz)
    

    sh2=subplot(224);
    sh2.Position(1) = sh2.Position(1) - 0.03;
    hold off
    scatter(d(:)./365,X(:),15,'k','filled','MarkerFaceAlpha',0.1,'MarkerEdgeColor','none');
    hold on
    scatter(d(:)./365,jk(:).*L./2,15,'r','filled','MarkerFaceAlpha',0.1,'MarkerEdgeColor','none');
    xlabel('Known divergence  (years)')
    ylabel('Estimated divergence (generations)')
    axis square
    box on
    title('(d) Clock vs. time')
    set(gca,'FontSize',fntsz)
    xlim([0 ceil(max(d(:)./365))])
    
    set(gcf,'Color','w')
    drawnow
end

if savefigs
    export_fig -r450 -opengl distances.png
    print -depsc -tiff -r300 -painters distances.eps
end
%%

if bingene_exists

    fig808=figure(808);
    fig808.Position = [152 350 350 350];
    clf
    
    % compare exact distances with estimates
    scatter(dgen(:),X(:),15,'k','filled','MarkerFaceAlpha',0.1,'MarkerEdgeColor','none');
    hold on
    scatter(dgen(:),jk(:).*L./2,15,'r','filled','MarkerFaceAlpha',0.1,'MarkerEdgeColor','none');
    
    xi = linspace(0,max(dgen(:)),100); % true number of generations
    xi_mut=xi.*eco_params.pneutral;    % true number of mutations
    psat = 1/2.*(1 - exp(-4.*(xi_mut./L))); % estimated saturated hamming distance
    
    ijk = linspace(0,i,100);
    xjk = psat;
    yjk = -1/2 * log(1 - xjk * 2/1);
    ejk = sqrt(xjk.*(1-xjk)./(L*((1-2.*xjk).^2)));
    X = p.*L./ eco_params.pneutral./2;
    pl1 =plot( xi,psat.*L./ eco_params.pneutral./2,'-','LineWidth',2,'Color',[0.5 0.5 0.5]);

    pl2 =plot(xi, yjk     .*L./2./ eco_params.pneutral,'-','LineWidth',2,'Color',[0.75 0 0]);
    pl2a=plot(xi,(yjk+ejk).*L./2./ eco_params.pneutral,':','LineWidth',2,'Color',[0.75 0 0]);
    pl2b=plot(xi,(yjk-ejk).*L./2./ eco_params.pneutral,':','LineWidth',2,'Color',[0.75 0 0]);
   
    xlim([0 max(dgen(:))]);
    xlabel('Known divergence (generations)')
    ylabel('Estimated divergence (generations)')
    axis square
    box on
    title('(c) Clock vs. generations')
    legend([pl2 pl1],'Corrected','Uncorrected','Location','SouthEast')
    set(gca,'FontSize',fntsz)
    

%     sh2=subplot(224);
%     sh2.Position(1) = sh2.Position(1) - 0.03;
%     hold off
%     scatter(d(:),X(:),15,'k','filled','MarkerFaceAlpha',0.1,'MarkerEdgeColor','none');
%     hold on
%     scatter(d(:),jk(:).*L./2,15,'r','filled','MarkerFaceAlpha',0.1,'MarkerEdgeColor','none');
%     xlabel('Known divergence  ( generations)')
%     ylabel('Estimated divergence (generations)')
%     axis square
%     box on
%     title('(d) Clock vs. time')
%     set(gca,'FontSize',fntsz)
%     xlim([0 ceil(max(d(:)))])
    
end

if savefigs
    set(gcf,'Color','w')
    drawnow
    export_fig -r450 -opengl clock_test.png
    print -depsc -tiff -r300 -painters clock_test.eps
end


return