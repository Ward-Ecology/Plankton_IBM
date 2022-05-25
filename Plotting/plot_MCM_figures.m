clear 
addpath Functions
bioinf_fcns  = bioinf_functions; 

outdir = 'output';
cd(outdir)
load('MCM_output.mat')

fntsz = 14;

nsample = 100; % max nsample terminal leaf nodes
nsample = min(nsample,size(biomass,1));

savefigs=true;


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
    exportgraphics(f999,'rates.png','Resolution',450)
end

return
