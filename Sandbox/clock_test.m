clear
clc

% add function path
addpath ../Functions
bioinf_fcns  = bioinf_functions; 

% option to save figures
savefigs = true;

npop   = 25; % population size 
ngenes = 50; % number of binary genes
pmut   = 1;  % probability of a single-point mutation at each generation

nsteps = ngenes.*53./pmut; % number of steps in random walk


ndists = nchoosek(npop+1,2)-npop; % number of paired distances between npop individuals
L      = ngenes*53;               % length of binary genome
tsteps = 1:nsteps;                % vector of timesteps

genes = uint64(zeros(npop,ngenes)); % array for discrtete binary genomes
rgb   =        zeros(npop,ngenes);  % array for continuous 'rgb' genomes
pbin  = zeros(nsteps,ndists);       % array for paired distances between binary genomes
prgb  = zeros(nsteps,ndists);       % array for paired distances between rgb genomes

f=waitbar(0,'Iterating...'); % open progress bar
for i=tsteps
    rgb       = bioinf_fcns.rgb_mutate(rgb,ones(size(rgb,1),1)); % mutate rgb gene
    prgb(i,:) = pdist(rgb); % calculate paired distances between rgb genomes
    
    if rand(1)<pmut % mutate with probability p
        genes     = bioinf_fcns.gene_mutate(genes,ones(size(genes,1),1)); % mutate binary genome
    end
    genestr   = bioinf_fcns.print_genomes(genes); % write binary genes out as ones and zeros
    pbin(i,:) = pdist(genestr-'0','Hamming'); % calculate paired distances between binary genomes

    waitbar(i./nsteps,f,'Iterating...'); % update progress bar
end
close(f) % close progress bar




%% Plot Figure A.1
f101 = figure(101);
f101.Position = [48 846 1091 491];
clf
fntsz=16;

%% RGB Genome
subplot(121) 
% plot accrued Euclidean distance in between rgb genes for 'ndists' pairs of individuals
plot(tsteps,prgb,'-','LineW',1,'Color',[0 0 0 0.01])
hold on

% convert number of steps to expected distance +/- 1 standard deviation
cfactor = gamma((ngenes+1)/2)./ gamma(ngenes/2); % statistical correction factor for n dimensions
d_expected = sqrt(4.*tsteps).*cfactor; % Equation A.7
d_err=sqrt(tsteps); % Equation A.8

% Plot expected distances between rgb genomes +/- 1 standard deviation
p1=plot([0 tsteps], [0 d_expected]     ,'-','LineWidth',2,'Color','k');    % plot expected distance 
plot([0 tsteps], [0 d_expected-d_err]     ,':','LineWidth',2,'Color','k'); % plot expected distance - 1 std dev
plot([0 tsteps], [0 d_expected+d_err]     ,':','LineWidth',2,'Color','k'); % plot expected distance + 1 std dev

% convert actual accrued Euclidean distances in between rgb genes to estimated number of steps based on theory
% i.e. accounting for non-linearity of expected distance as function of time/generations
t_est = prgb.^2./(4.*cfactor.^2); % solve equation A.7 for number of generations
plot(tsteps,t_est,'-','LineW',1,'Color',[1 0 0 0.01])

% convert expected accrued distances (+/- 1 s.d.) to expected estimates +/- error 
y=d_expected.^2./(4.*cfactor.^2);            % 1:1 relationship
ym = (d_expected-d_err).^2./(4.*cfactor.^2); % 1:1 relationship - expected error
yp = (d_expected+d_err).^2./(4.*cfactor.^2); % 1:1 relationship + expected error

p2=plot([0 tsteps], [0 y],'-','LineWidth',2,'Color','r'); % plot 1:1 relationship
plot([0 tsteps], [0 ym],':','LineWidth',2,'Color','r');   % plot 1:1 relationship - expected error
plot([0 tsteps], [0 yp],':','LineWidth',2,'Color','r');   % plot 1:1 relationship + expected error

% format subplot
axis([0 nsteps 0 nsteps.*1.125])
set(gca,'XTick',0:800:nsteps*1.125,'YTick',0:800:nsteps*1.125)
title('(a) RGB genome')
xlabel('True number of generations')
ylabel('Distance or estimated number of generations')
legend([p1 p2],'p (distance)','estimated t_{gen}','Location','NorthWest')
set(gca,'FontSize',fntsz)

%% Binary Genome

subplot(122)

% plot accrued Euclidean distance in between binary genes for 'ndists' pairs of individuals
plot(tsteps,pbin.*L,'-','LineW',1,'Color',[0 0 0 0.01])
hold on

% convert number of steps to expected distance +/- 1 standard deviation
dbin = 1/2.*(1 - exp(-4.*tsteps.*pmut./L)); % estimated saturated  hamming distance (normalised)
derr = sqrt(dbin.*(1-dbin)./L);

% Plot expected distances between binary genomes +/- 1 standard deviation
p1 =plot(tsteps, dbin.*L     ,'-','LineWidth',2,'Color','k');
plot(tsteps, (dbin-derr).*L ,':','LineWidth',2,'Color','k');
plot(tsteps, (dbin+derr).*L ,':','LineWidth',2,'Color','k');

% convert actual accrued Hamming distances in between binary genes to estimated number of steps based on theory
% Using 2-base Jukes-Cantor
jk = -1/4 * log(1 - 2.*pbin) ./pmut;
plot(tsteps,jk.*L,'-','LineW',1,'Color',[1 0 0 0.01])

% convert expected accrued distances (+/- 1 s.d.) to expected estimates +/- error 
yjk = -1/4 * log(1 - 2.*dbin)./pmut;
ejk = sqrt(dbin.*(1-dbin)./(4*L*((1-2.*dbin).^2)))./pmut;

% convert expected accrued distances (+/- 1 s.d.) to expected estimates +/- error 
p2=plot(tsteps,yjk.*L,'-','LineWidth',2,'Color','r');
plot(tsteps,(yjk+ejk).*L,':','LineWidth',2,'Color','r');
plot(tsteps,(yjk-ejk).*L,':','LineWidth',2,'Color','r');
  
% format subplot
axis([0 nsteps 0 nsteps.*1.125])
set(gca,'XTick',0:800:nsteps*1.125,'YTick',0:800:nsteps*1.125)
title('(b) Binary genome')
xlabel('True number of generations')
ylabel('Distance or estimated number of generations')
legend([p1 p2],'p (distance)','estimated t_{gen}','Location','NorthWest')
set(gca,'FontSize',fntsz)

%%
if savefigs
    set(gcf,'Color','w')
    exportgraphics(f101,'clock_test.png','Resolution',450)
end



























