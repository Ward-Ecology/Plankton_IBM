clear 
addpath ~/GitHub/Plankton_IBM/Functions
bioinf_fcns  = bioinf_functions; 

outdir = '~/GitHub/Plankton_IBM/output/';
cd(outdir)
load('MCM_output.mat')

fntsz = 14;

nsample = 100; % max nsample terminal leaf nodes
nsample = min(nsample,size(biomass,1));

savefigs=true;

%% sample MCM
biomass(biomass<0)=0;
ind_MCM = sort(randsample(numel(biomass),nsample,true,biomass));
gensample = genome(ind_MCM,:);
% convert condensed format to full binary strings
genome_string = bioinf_fcns.print_genomes(gensample);

% generate hamming-distance matrix
p = pdist(genome_string-'0','Hamming');
p = squareform(p); % make square matrix
L = eco_params.ngenome.*64;
X = p.*L./ eco_params.pneutral./2; % convert from fraction of bits different to total bits different
% also account for mutation probability
% 2-base jukes-cantor correction
pcorr=p;
pcorr(pcorr>0.5)=0.5;
jk = -1/2 * log(1 - pcorr * 2/1);

%  normalise rgb gene (centered on origin)
plotclr = rgb(ind_MCM,:);
maxclr  = max(abs(plotclr(:)));
minclr  = -maxclr;
plotclr = (plotclr-minclr)./(maxclr-minclr);

%%
f1=figure(1);
f1.Position=[1000 561 772 777];
clf

% plot biomass contours (phenotype through time)

subplot(311)
bio_out=yout(:,2:end-1)';
perc=100.*bio_out./sum(bio_out,1);
ycoord = linspace(eco_params.T_bins(1),eco_params.T_bins(end),size(perc,1));

contourf(tout./env_forcing.daysperyear,ycoord,perc,0:10,'LineColor','none');
axis xy
cmap = parula(10);
cmap(1,:)=[1 1 1];
colormap(cmap)
caxis([0 10])
ch=colorbar('Location','West');
ch.Position(4)=ch.Position(4).*0.4;
ch.Ticks=0:5:10;
ch.FontSize=fntsz;
hold on
% plot(tout./env_forcing.daysperyear,env_forcing.Tfunc(tout),'color',[0.5 0.5 0.5],'LineWidth',1)
title('(a) Population dynamics','FontSize',fntsz)
xlabel('Time (years)','FontSize',fntsz)
ylabel('Temperature (°C)','FontSize',fntsz)

% plot rgb colours on same axes
subplot(312)

rgb_perm = permute(rgb_out,[2 1 3]);
lo_prc = prctile(rgb_perm(:),5);
hi_prc = prctile(rgb_perm(:),95);
rgb_norm = (rgb_perm-lo_prc)./(hi_prc-lo_prc);
im=image(tout./env_forcing.daysperyear,ycoord,rgb_norm(:,:,1:3));
im.AlphaData = bio_out./max(bio_out,[],1);
hold on
contour(tout./env_forcing.daysperyear,ycoord,perc,[0 1],'LineColor','k');
% plot(tout./env_forcing.daysperyear,env_forcing.Tfunc(tout),'color',[0.5 0.5 0.5],'LineWidth',1)
axis xy
box on
title('(b) rgb genomes (opacity scaled by biomass)','FontSize',fntsz)
xlabel('Time (years)','FontSize',fntsz)
ylabel('Temperature (°C)','FontSize',fntsz)

subplot(313)
im=image(tout./env_forcing.daysperyear,ycoord,rgb_norm(:,:,1:3));axis xy
hold on

% contour(tout./env_forcing.daysperyear,ycoord,perc,[0 1],'LineColor','k');
% contour(tout./env_forcing.daysperyear,ycoord,mu_out-eco_params.gamma_0_mort,[0 0],'LineColor','k');
plot(tout./env_forcing.daysperyear,env_forcing.Tfunc(tout),'color',[0.5 0.5 0.5],'LineWidth',1)

box on
title('(c) rgb genomes (unscaled)','FontSize',fntsz)
xlabel('Time (years)','FontSize',fntsz)
ylabel('Temperature (°C)','FontSize',fntsz)


set(gcf,'Color','w')



if savefigs
    export_fig -r450 time_series.png
    print -depsc -tiff -r300 -painters time_series.eps
end
%%
clear bingenes

for i=1:nsample
    bingenes(i).Sequence=char(97+genome_string(i,:)-'0');
    bingenes(i).Header='';
end
 
% distances = seqpdist(bingenes,'Method','Jukes-Cantor','Alpha','DNA').*L;
UPGMAtree = seqlinkage(jk.*L,'UPGMA',bingenes);
NJtree    = seqneighjoin(jk.*L,'equivar',bingenes);

h1 = plot(UPGMAtree,'orient','left');
% axis xy
title('UPGMA Distance Tree of agents using Jukes-Cantor model');
ylabel('Evolutionary distance')
if savefigs
    print -depsc -tiff -r300 -painters UPGMA.eps
end

h2 = plot(NJtree,'orient','left');
% axis xy
title('Neighbor-Joining Distance Tree of agents using Jukes-Cantor model');
ylabel('Evolutionary distance')
if savefigs
    print -depsc -tiff -r300 -painters UPGMA.eps
end




%Comparing Tree Topologies
% Notice that different phylogenetic reconstruction methods result in
% different tree topologies. The neighbor-joining tree groups Chimp
% Vellerosus in a clade with the gorillas, whereas the UPGMA tree groups it
% near chimps and orangutans. The |getcanonical| function can be used to
% compare these isomorphic trees.

sametree = isequal(getcanonical(UPGMAtree), getcanonical(NJtree));


XData=[h1.LeafDots.XData h1.BranchDots.XData];
YData=[h1.LeafDots.YData h1.BranchDots.YData];
% leafOrder = optimalleaforder(Z,jk);

[adj,id,~]=getmatrix(UPGMAtree);

leafnames   = strvcat(get(UPGMAtree,'LeafNames'));
sortnodes = str2num(leafnames(:,5:end));


%%


f2=figure(2);
f2.Position = [139 103 1239 468];
clf
% % Create hierarchical cluster tree
% Z = linkage(jk.*L./2,'average');
% leafOrder = optimalleaforder(Z,jk);

% plot dendrogram
sh1=subplot(161);
G = digraph(adj);

% add extra nodes to square up dendrogram
XDta = XData;
YDta = YData;
c = size(G.Nodes,1); 
for i=1:size(G.Nodes,1);
    preID = predecessors(G,i);
    if numel(preID)>0
        c =c+1;
        G = addedge(G,preID,c);
        XDta=[XDta XDta(preID)];
        YDta=[YDta YDta(i)];
        
        G = addedge(G,c,i);
        G = rmedge(G,preID,i);
    end
end
        
    

g=plot(G);

g.XData=XDta;
g.YData=YDta;
xlim
% xlim([0 450])

g.LineWidth=2;
g.NodeColor='none';
g.ArrowSize=0;
g.EdgeAlpha=1;
g.NodeLabel=[];
g.EdgeColor='k';


% [~,ind]=sort(outperm);
hold on
axis off
ylim([0.5 nsample+0.5])

genmat=genome_string-'0';
rgbmat=zeros(nsample,3200,3);
for i=1:nsample
    clr=plotclr(i,:)';
    rgbmat(i,:,:) = (clr*genmat(i,:))';
end
% rgbmat(rgbmat==0)=1;


subplot(1,6,[2 6])
imagesc(rgbmat(sortnodes,:,:))
axis xy
axis tight
box on
set(gca,'XTick',[],'YTick',[])
xlim([1 128])

cmap=redblue(numel(biomass));
ch=colorbar('Location','eastoutside','Ticks',[],'FontSize',fntsz);
ylabel(ch, 'Thermal optimum')
colormap(cmap(ind_MCM(sortnodes),:))

sh1.Position=[0.15  0.1100    0.1023    0.8150];

set(gcf,'Color','w')


if savefigs
    export_fig -r450 genome.png
    print -depsc -tiff -r300 -painters genome.eps
end
%%

dmat=jk.*L./2;
titl='Estimated divergence';

fig3=figure(3);
fig3.Position = [152 700 500 500];
clf

imagesc(dmat(sortnodes,sortnodes));
axis square
axis xy
ch=colorbar('FontSize',fntsz);
ylabel(ch, 'Divergence')
ax1_pos = get(gca,'Position'); % position of first axes
set(gca,'Position',ax1_pos,'XTick',[],'YTick',[]);
caxis([0 max(dmat(:))])
caxis([0 450]);
title(titl,'FontSize',fntsz)
colormap(parula(64));

% rgb gene
ax2 = axes('Position',ax1_pos,'XTick',[],'YTick',[]);
axis square
axis off
ch=colorbar('Location','southoutside','Ticks',[],'FontSize',fntsz);
ylabel(ch, 'RGB gene')
colormap(ax2,plotclr(sortnodes,:))
ax2.Position = ax1_pos;

% thermal optima
cmap=redblue(numel(biomass));
ax3 = axes('Position',ax1_pos,'XTick',[],'YTick',[]);
axis square
axis off
ch=colorbar('Location','westoutside','Ticks',[],'FontSize',fntsz);
ylabel(ch, 'Thermal optimum')
colormap(ax3,cmap(ind_MCM(sortnodes),:))
ax3.Position = ax1_pos;

set(gcf,'Color','w')


if savefigs
    export_fig -r450 distances.png
    print -depsc -tiff -r300 -painters distances.eps
end

%%


% get environmental forcing (i.e. temperature)
T  = env_forcing.Tfunc(env_forcing.tmax);
% get state variables from V_vect
N  = yout(end,1);
P  = yout(end,2:eco_params.jpmax+1)';
D  = yout(end,eco_params.jpmax+2:end);
% evaluate limitation factors
F_T = Temperature_function(T,eco_params)';
F_N = Nutrient_function(N,eco_params);
% remame parameters in shorter form
mu0   = eco_params.mu0_max_rate;
gamma = eco_params.gamma_0_mort;
tau0  = eco_params.tau_remin;
kN    = eco_params.kN_half_sat;
% ecosystem model equations
%==========================
% ODEs
mu  = mu0.*F_T.*F_N;
muP  = mu.*P;
Tmut = eco_params.TT.*mu';
dPdt = + muP + Tmut*P - gamma .* P;


Nstar = gamma.*kN./(mu0-gamma);
T1000 = linspace(13,17,1000);
munet = mu0.*exp( -((15 - T1000)./6).^2 ).*Nutrient_function(Nstar,eco_params) - gamma;


f999 = figure(999);
f999.Position = [165   651   560   651];
clf
subplot(3,4,[1 3])
plot(eco_params.T_optimum,P./eco_params.b_0,'k','LineW',2)
ylabel('Abundance (individuals m^{-3})');
set(gca,'XTick',13:17);
xlim([13 17])
title('(a) Equilibrium trait distribution')

subplot(3,4,[5 7])
[ax,h1,h2] = plotyy(T1000,munet,T1000,-1./munet,@plot,@semilogy);
ax(1).YLim=[-0.01 0];
ax(1).YTick=-0.01:0.005:0;
ax(1).XLim=[13 17];
ax(1).YColor='k';
h1.Color='k';
h1.LineWidth = 2;
h1.LineStyle = '--';
ylabel(ax(1), 'Fitness (d^{-1})');
ax(1).XTick=13:17;
set(ax(1),'Box','off')
title('(b) Fitness and exclusion time scales')

ax(2).YLim=[30 360000];
ax(2).YTick=[30 360.*[1 10 100 1000]];
ax(2).YTickLabel = {'Month','Year','Decade','Century','Millenium'}	;
ax(2).YMinorTick = 'off';
ax(2).YGrid = 'on';
ax(2).YMinorGrid = 'off';
ax(2).XLim=[13 17];
ax(2).YColor='k';
h2.Color='k';
h2.LineWidth = 2;
ylabel(ax(2), 'Exclusion time scale');

subplot(3,4,[9 11])
plot(eco_params.T_optimum,dPdt./eco_params.b_0,'k','LineW',2,'LineStyle','-')
hold on
plot(eco_params.T_optimum,(muP - gamma.*P)./eco_params.b_0,'k','LineW',2,'LineStyle','--')
plot(eco_params.T_optimum,(Tmut*P)./eco_params.b_0,'k','LineW',2,'LineStyle',':')
legend('net growth','birth-death','mutation')
xlim([13 17])
ylabel('Overall rates (Individuals m^{-3} day^{-1})')
xlabel('Thermal optima (^\circC)')
set(gca,'XTick',13:17);
title('(c) Biomass sources and sinks')

drawnow

set(gcf,'Color','w')
if savefigs
    export_fig -r450 rates.png
    print -depsc -tiff -r300 -painters rates.eps
end

return
%% Animated bar chart 

for dy=1:numel(tout)
    barclr = squeeze(rgb_norm(:,dy,:));
    
%     b = bar(ycoord,100.*bio_out(:,dy)./sum(bio_out(:,dy)));
%     b = bar(ycoord,log10(100.*bio_out(:,dy)./sum(bio_out(:,dy))));
    hold off
    b = bar(ycoord,10.*ones(size(bio_out(:,dy))));
    hold on
    plot(ycoord,100.*bio_out(:,dy)./sum(bio_out(:,dy)),'k','LineW',2)
    
    b.BarWidth = 1;
    b.FaceColor = 'flat';
    b.CData = barclr;
    axis([min(ycoord) max(ycoord) 0 10])
    xlabel('Thermal optimum (^\circC)')
    ylabel('Biomass fraction (%)')
    title(num2str(dy./env_forcing.daysperyear))
    drawnow
end


