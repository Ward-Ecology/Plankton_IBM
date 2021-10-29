clear 
addpath ~/GitHub/Plankton_IBM/Functions
bioinf_fcns  = bioinf_functions; 

outdir = '~/GitHub/Plankton_IBM/output/';
cd(outdir)

fntsz = 14;

savefigs=true;


IBM_output=load('IBM_output.mat');

eco_params  = IBM_output.eco_params;
env_forcing = IBM_output.env_forcing;
if eco_params.use_bioinf
    genome      = IBM_output.bioinformatics.genome;
    rgb         = IBM_output.bioinformatics.rgb;
end

t_opt       = eco_params.T_optimum;

%%
bio_out=IBM_output.yout(:,2:end-1);
perc=(100.*bio_out./sum(bio_out,2));
ycoord = linspace(eco_params.T_bins(1),eco_params.T_bins(end),size(perc,2));

ylimits = [8 22];

f102=figure(102);
f102.Position=[91 900 600 400];

clf
% imagesc(IBM_output.tout,ycoord,perc');
cmap=flipud(gray(8));
contourf(IBM_output.tout,ycoord,perc',[0:8],'LineColor',[1 1 1].*0.2,'LineWidth',0.25);
hold on
caxis([0 8])
colormap(cmap)

hold on
nyears=env_forcing.tmax./env_forcing.daysperyear;
axis([0 env_forcing.tmax ylimits])
set(gca,'XTick',linspace(0,env_forcing.tmax,min(nyears+1,21)));
set(gca,'XTickLabel',linspace(0,nyears,min(nyears+1,21)));
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


T=0:eco_params.delta_t:env_forcing.tmax;
scatter(T,env_forcing.Tfunc(T),1,'r','filled','MarkerEdgeAlpha',0.1)


if savefigs
    export_fig -r450 time_series.png
    print -depsc -tiff -r300 -painters time_series.eps
end
%%

%%

% N = eco_params.gamma_0_mort.*eco_params.kN_half_sat./(eco_params.mu0-eco_params.gamma_0_mort);
N = mean(IBM_output.yout(end-359:end,1)./eco_params.V,1);
% P = mean(bio_out(end-359:end,:)'./eco_params.V,2);


T   = env_forcing.Tfunc(1:env_forcing.tmax);
T_phen = ycoord';

F_T = exp( -((T - T_phen)./eco_params.theta_norm).^2 );
F_N = Nutrient_function(N,eco_params)';

% % Transport matrix
% eco_params.deltaM = eco_params.sigmaM.^2 ...      % genotype transfer function
%             ./ ( 3.*eco_params.deltaT_opt.^2 );
% n      = eco_params.jpmax;                      % n 'genotypes'
% delta  = ones(n,1).*eco_params.delta;          % expanded vector of genotype transfer function
% TT     = spdiags([delta delta],[-1 1], n, n); % set off-diagonals to delta M
% TT     = spdiags(-sum(TT,2),0, TT);             % set diagonal to conserve mass
% TT = TT + speye(size(TT));

fitness =    (eco_params.mu0 .* F_T .* F_N) - eco_params.gamma_0_mort;
% fitness = TT*(eco_params.mu0 .* F_T .* F_N) - eco_params.gamma_0_mort;

f654 = figure(654);
f654.Position=[91 900 600 400];

clf
% sh1=subplot(211);
imagesc(IBM_output.tout,T_phen,fitness)
hold on
contour(IBM_output.tout,ycoord,fitness,-[1 1]./(15.*360),'LineColor','k','LineWidth',1);
colormap(flipud(redblue))
ch=colorbar;
caxis([-0.1 0.1])
set(gca,'XTick',linspace(0,env_forcing.tmax,nyears+1));
set(gca,'XTickLabel',linspace(0,nyears,nyears+1));
% xlim(env_forcing.tmax-[15.*360 0])
title('Invasion fitness (1/day)')
% ch.TickLabels = cellstr(num2str(1./abs(ch.Ticks'),'%3.0f'));
% ch.TickLabels(ticks>360) = cellstr(num2str(ticks(ticks>360)'./360));

ylim(ylimits)
% 
% sh2=subplot(212);
% ticks = [1 7 30 [1 2 5 10].*360];
% excl_time = 1./abs(fitness);
% excl_time(fitness>=0)=max(ticks);
% imagesc(IBM_output.tout,T_phen,log10(excl_time))
% hold on
% contour(IBM_output.tout,ycoord,perc',[1 1],'LineColor',[1 1 1].*0.2,'LineWidth',0.25);
% colormap(sh2,parula)
% ch=colorbar;
% ch.Ticks = log10(ticks);
% ch.TickLabels = cellstr(num2str(ticks'));
% ch.TickLabels(ticks>360) = cellstr(num2str(ticks(ticks>360)'./360));
% caxis(log10([ticks(1) ticks(end)]))
% set(gca,'XTick',linspace(0,env_forcing.tmax,nyears+1));
% set(gca,'XTickLabel',linspace(0,nyears,nyears+1));
% xlim(env_forcing.tmax-[360 0])
% title('Exclusion timescale (days)')


set(gcf,'Color','w')
if savefigs
    export_fig -png -r450 -opengl fitness.png
    print -depsc -r300 -painters fitness.eps
end


%% sample terminal leaf nodes 
nsample = 1000; % max nsample terminal leaf nodes
iextant = find(IBM_output.y(2:end-1));
nextant = numel(iextant);
nsample = min(nsample,nextant); % decrease if not enough nodes
isample = randsample(iextant,nsample); % sample without replacement

[~,ii]=sort(t_opt(isample));

isample = isample(ii);

%% bioinformatics toolbox
genome_string = bioinf_fcns.print_genomes(genome(isample,:));

 
% generate hamming-distance matrix
L = eco_params.ngenome.*53;
p = pdist(genome_string-'0','Hamming');
p = squareform(p); % make square matrix
X = p.*L./ eco_params.pneutral./2; % convert from fraction of bits different to total bits different
% also account for mutation probability
% 2-base jukes-cantor correction
pcorr=p;
pcorr(pcorr>0.5)=NaN;
jk = -1/2 * log(1 - pcorr * 2/1);
jk = jk./eco_params.pneutral;
    
figure(333)
hist(p(:),100)
    
clear bingenes 
genmat=genome_string-'0';
for i=1:nsample
    bingenes(i).Sequence=char(97+genmat(i,:)-'0');
    bingenes(i).Header=num2str(i);
end


UPGMAtree = seqlinkage(jk.*L,'UPGMA',bingenes);

leafnames   = strvcat(get(UPGMAtree,'LeafNames'));
sortnodes = str2num(leafnames);
[~,isort] = sort(sortnodes);
UPGMAtree = reorder(UPGMAtree,isort,'approximate',true);

nclusters=20;
[clus,nclus,steps] = cluster(UPGMAtree,[],'maxclust',nclusters);


h1 = plot(UPGMAtree,'orient','left');
for i=1:nclusters
   set(h1.BranchLines(nclus==i),'Color',rand(1,3))
end
% axis xy
title('UPGMA Distance Tree of agents using Jukes-Cantor model');
ylabel('Evolutionary distance')
if savefigs
    print -depsc -tiff -r300 -painters UPGMA.eps
end


Xdist=[h1.LeafDots.XData h1.BranchDots.XData];
Ydist=[h1.LeafDots.YData h1.BranchDots.YData];
bclr = reshape([h1.BranchLines.Color],3,[])';
% Xdist=[h2.LeafDots.XData h2.BranchDots.XData];
% Ydist=[h2.LeafDots.YData h2.BranchDots.YData];

%%
f303=figure(303);
f303.Position = [139 103 1200 400];
clf

tree = UPGMAtree;

sh2=subplot(161);
adj = getmatrix(tree);
g2 = digraph(adj);
% add extra nodes to square up forks in dendrogram
XDta = Xdist;
YDta = Ydist;
c = size(g2.Nodes,1); 
nclus2 = nclus';
for i=1:size(g2.Nodes,1);
    preID = predecessors(g2,i);
    if numel(preID)>0
        c =c+1;
        g2 = addedge(g2,preID,c);
        XDta=[XDta XDta(preID)];
        YDta=[YDta YDta(i)];
        nclus2=[nclus2 nclus2(i)];
            
        g2 = addedge(g2,c,i);
        g2 = rmedge(g2,preID,i);
    end
end
g=plot(g2);
g.XData=XDta;
g.YData=YDta;

g2.Edges.EndNodes(:,2);

% xlim([0 450])

g.LineWidth=1;
g.NodeColor='none';
g.ArrowSize=0;
g.EdgeAlpha=1;
g.NodeLabel=[];

cmap=rand(nclusters,3);
g.EdgeColor=cmap(nclus2(g2.Edges.EndNodes(:,2)),:);
axis ij
% [~,ind]=sort(outperm);
hold on
axis off
ylim([0.5 nsample+0.5])
title({'(a) phylogeny'},'FontSize',fntsz)


%  normalise rgb gene (centered on origin)
rgb_samp=rgb(isample,:);

maxclr  = max(abs(rgb_samp));
minclr  = -maxclr;
plotclr = (rgb_samp-minclr)./(maxclr-minclr);

leafnames   = strvcat(get(tree,'LeafNames'));
sortnodes = flipud(str2num(leafnames));
% Group thermal optima into bins
cbin = discretize(t_opt(isample(sortnodes)),linspace(min(t_opt(iextant)),max(t_opt(iextant)),64+1));

genmat=genome_string(sortnodes,:)-'0';
rgbmat=zeros(nsample,size(genmat,2),3);
for i=1:nsample
    clr=plotclr(sortnodes(i),:)';
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
xlim([1 2^8])

cmap=redblue(64);
ch=colorbar('Location','eastoutside','Ticks',[],'FontSize',fntsz);
ylabel(ch, 'Thermal optimum','FontSize',fntsz)
colormap(cmap(cbin,:))
title('(b) neutral genome','FontSize',fntsz)

sh2.Position=[0.15  0.1100    0.1023    0.8150];
set(gcf,'Color','w')

if savefigs
    export_fig -r450 genome.png
    print -depsc -tiff -r300 -painters genome.eps
end

%%




figure(111)
hold on
xx=sort(Xdist);
xx=xx-max(xx);
xx=xx(1:nsample);

plot(xx,1:nsample,'LineWidth',2)

thalf = 2.*nextant/nsample;
plot(xx,nsample./(1-xx./thalf),'k','LineWidth',2);
xlim([-10.^floor(log10(nextant)) -1])

xlabel('Generations before end of simulation')
ylabel('Coalescences')

if savefigs
    export_fig -r450 coalesence.png
    print -depsc -tiff -r300 -painters coalesence.eps
end
%%

dmat=jk(sortnodes,sortnodes).*L./2;
titl='Estimated divergence';

fig1=figure(202);
fig1.Position = [152 700 400 400];
clf
cmap=redblue(64);

imagesc(dmat);
axis square
axis xy
ch=colorbar('FontSize',fntsz);
ylabel(ch, 'Divergence (generations)')
ax1_pos = get(gca,'Position'); % position of first axes
set(gca,'Position',ax1_pos,'XTick',[],'YTick',[]);
caxis([0 max(dmat(:))])
cx=caxis;
title(titl,'FontSize',fntsz)

% rgb gene
ax2 = axes('Position',ax1_pos,'XTick',[],'YTick',[]);
axis square
axis off
ch=colorbar('Location','southoutside','Ticks',[],'FontSize',fntsz);
ylabel(ch, 'RGB gene')
colormap(ax2,plotclr(sortnodes,:))
ax2.Position = ax1_pos;

% thermal optima

ax3 = axes('Position',ax1_pos,'XTick',[],'YTick',[]);
axis square
axis off
ch=colorbar('Location','westoutside','Ticks',[],'FontSize',fntsz);
ylabel(ch, 'Thermal optimum')
colormap(ax3,cmap(cbin,:))
ax3.Position = ax1_pos;

set(gcf,'Color','w')
if savefigs
    export_fig -r450 distances.png
    print -depsc -tiff -r300 -painters distances.eps
end
