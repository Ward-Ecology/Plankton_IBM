clear
clc

addpath Functions
bioinf_fcns  = bioinf_functions; 

savefigs = true;

npop   = 25;
ngenes = 100;
pmut   = 1;

% number of steps in random walk
nsteps = ngenes.*53./pmut;

ndists = nchoosek(npop+1,2)-npop;
L      = ngenes*53;
tsteps = 1:nsteps;

genes = uint64(zeros(npop,ngenes));
rgb   =        zeros(npop,ngenes);
pbin  = zeros(nsteps,ndists);
prgb  = zeros(nsteps,ndists);

cfactor       = gamma((ngenes+1)/2)./ gamma(ngenes/2); % correction factor for n dimensions

clf
for i=tsteps
    rgb       = bioinf_fcns.rgb_mutate(rgb,ones(size(rgb,1),1));
    prgb(i,:) = pdist(rgb);
    
    if rand(1)<pmut
        genes     = bioinf_fcns.gene_mutate(genes,ones(size(genes,1),1));
    end
    genestr   = bioinf_fcns.print_genomes(genes);
    pbin(i,:) = pdist(genestr-'0','Hamming');
    clc,disp(i)
end


%%
clf
fntsz=16;

subplot(121)
% plot accrued Euclidean distance through time
plot(tsteps,prgb,'-','LineW',1,'Color',[0 0 0 0.01])
hold on
% plot(tsteps,mean(prgb,2),'g-','LineW',2) % plot evaluated mean
% plot(tsteps,mean(prgb,2)*[1 1]+std(prgb,[],2)*[-1 1] ,'g--','LineW',2) % plot +/- evaluated stdev

% convert number of steps to expected distance
d_expected = sqrt(4.*tsteps).*cfactor;
p1=plot([0 tsteps], [0 d_expected]     ,'-','LineWidth',2,'Color','k');
d_err=sqrt(tsteps);
plot([0 tsteps], [0 d_expected-d_err]     ,':','LineWidth',2,'Color','k');
plot([0 tsteps], [0 d_expected+d_err]     ,':','LineWidth',2,'Color','k');

% convert accrued distance to estimate number of steps
t_est = prgb.^2./(4.*cfactor.^2);
plot(tsteps,t_est,'-','LineW',1,'Color',[1 0 0 0.01])
% plot(tsteps,mean(t_est,2),'k-','LineW',2) % plot evaluated mean
% plot(tsteps,mean(t_est,2)*[1 1]+std(t_est,[],2)*[-1 1] ,'k--','LineW',2) % plot +/- evaluated stdev


y=d_expected.^2./(4.*cfactor.^2);
p2=plot([0 tsteps], [0 y],'-','LineWidth',2,'Color','r');
yp = (d_expected+d_err).^2./(4.*cfactor.^2);
ym = (d_expected-d_err).^2./(4.*cfactor.^2);

plot([0 tsteps], [0 ym],':','LineWidth',2,'Color','r');
plot([0 tsteps], [0 yp],':','LineWidth',2,'Color','r');

axis([0 nsteps 0 nsteps.*1.125])
set(gca,'XTick',0:800:nsteps*1.125,'YTick',0:800:nsteps*1.125)
title('(a) RGB genome')
xlabel('True number of generations')
ylabel('Distance or estimated number of generations')
legend([p1 p2],'p (distance)','estimated t_{gen}','Location','NorthWest')
set(gca,'FontSize',fntsz)

%%

subplot(122)
cla
% plot accrued base differences through time
plot(tsteps,pbin.*L,'-','LineW',1,'Color',[0 0 0 0.01])
hold on
% plot(tsteps,mean(pbin,2).*L./2,'g-','LineW',1) % plot evaluated mean
% plot(tsteps,(mean(pbin,2)*[1 1]+std(pbin,[],2)*[-1 1]).*L./2 ,'g--','LineW',1) % plot +/- evaluated stdev

jk = -1/4 * log(1 - 2.*pbin) ./pmut;
% jke = sqrt(pbin.*(1-pbin)./(4*L*((1-2.*pbin).^2)));
plot(tsteps,jk.*L,'-','LineW',1,'Color',[1 0 0 0.01])

% estimate distances from generations
tgen = linspace(0,nsteps,100); % true number of generations
dbin = 1/2.*(1 - exp(-4.*tgen.*pmut./L)); % estimated saturated  hamming distance (normalised)
derr = sqrt(dbin.*(1-dbin)./L);

p1 =plot(tgen, dbin.*L     ,'-','LineWidth',2,'Color','k');
plot(tgen, (dbin-derr).*L ,':','LineWidth',2,'Color','k');
plot(tgen, (dbin+derr).*L ,':','LineWidth',2,'Color','k');


yjk = -1/4 * log(1 - 2.*dbin)./pmut;
ejk = sqrt(dbin.*(1-dbin)./(4*L*((1-2.*dbin).^2)))./pmut;


p2=plot( tgen,yjk.*L,'-','LineWidth',2,'Color','r');%%
plot(tgen,(yjk+ejk).*L,':','LineWidth',2,'Color','r');
plot(tgen,(yjk-ejk).*L,':','LineWidth',2,'Color','r');
  

  
axis([0 nsteps 0 nsteps.*1.125])
set(gca,'XTick',0:800:nsteps*1.125,'YTick',0:800:nsteps*1.125)
title('(b) Binary genome')
xlabel('True number of generations')
ylabel('Distance or estimated number of generations')
legend([p1 p2],'p (distance)','estimated t_{gen}','Location','NorthWest')
set(gca,'FontSize',fntsz)

if savefigs
    set(gcf,'Color','w')
    export_fig -r450 -opengl /Users/baw103/GitHub/RGB_evolution/Figures/clock_test.png
    print -depsc -tiff -r300 -painters /Users/baw103/GitHub/RGB_evolution/Figures/clock_test.eps
end

return
  %%
    
    
    jk = -1/2 * log(1 - p * 2/1);
    jke = sqrt(p.*(1-p)./(L*((1-2.*p).^2)));
    
    subplot(211)
    plot(i, jk.*L./2,'r.')
    hold on
    plot(i, p.*L./2,'k.')
    
 
    
    if exist('pl1')
        delete(pl1)
%         delete(pl1a)
%         delete(pl1b)
%         delete(pl2)
%         delete(pl2a)
%         delete(pl2b)
    end
    subplot(211)

    pl2 =plot(xi, yjk     .*L./2,'-','LineWidth',2,'Color',[0.75 0 0]);
    pl2a=plot(xi,(yjk+ejk).*L./2,':','LineWidth',2,'Color',[0.75 0 0]);
    pl2b=plot(xi,(yjk-ejk).*L./2,':','LineWidth',2,'Color',[0.75 0 0]);
    
%     axis([0 i 0 max(jk*L/2)])
    
    if rem(i,30)==0
        drawnow
    end



%%
clear

% expected RMS displacement of lattice or Gaussian random walk
walk_type = 'Gaussian'; % 'Gaussian' or 'Lattice'
% random walk in ndim dimensions
ndim  = 50;
% 
npop=2;
 % number of repetitions
niter = 1e2;


cfactor = gamma((ndim+1)/2)./ gamma(ndim/2);

d = zeros(nsteps,niter); % initialise distance array
for k=1:niter % repeat niter times
    coord=zeros(npop,ndim); % initialise walk at origin
    for t=1:nsteps % walk for nsteps
        % take random walk step
        switch walk_type
            case 'Lattice'
                coord = coord + sign(randn(npop,ndim));
            case 'Gaussian'
                coord = coord + randn(npop,ndim)
        end
        % save euclidean distance between two points at time t
        d(t,k)=pdist([coord],'euclidean');
    end
end
 
figure(1)
clf
title(['Euclidean distance after n timesteps in ' num2str(ndim) ' dimensions'])
tsteps=1:nsteps;
stairs(tsteps,d,'-','LineW',0.1,'Color',[0.75 0.75 0.75])
hold on
plot(tsteps,mean(d,2),'k-','LineW',2) % plot evaluated mean
plot(tsteps,mean(d,2)*[1 1]+std(d,[],2)*[-1 1] ,'k--','LineW',2) % plot +/- evaluated stdev
% calculate and plot theoretical mean (and stdev?)
d_expected = sqrt(2.*npop.*tsteps).*cfactor;
% d_stdev = ???
plot(tsteps,d_expected,'g--','LineW',2)


plot(tsteps,d_expected+sqrt(2.*npop.*tsteps.*0.5^2),'g:','LineW',2)
plot(tsteps,d_expected-sqrt(2.*npop.*tsteps.*0.5^2),'g:','LineW',2)


























