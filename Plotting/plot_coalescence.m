clear
addpath ~/GitHub/Plankton_IBM/Functions
bioinf_fcns  = bioinf_functions;

f1=figure(111);
clf

rundirs = {...
           'constant_mut10_IBM_1e4',...
           'constant_mut10_IBM_1e5',...
           'constant_mut10_IBM_1e6',...
           'constant_mut10_IBM_1e7',...
            ...
           'stepfunction1_mut10_IBM_1e4',...
           'stepfunction1_mut10_IBM_1e5',...
           'stepfunction1_mut10_IBM_1e6',...
           'stepfunction1_mut10_IBM_1e7',...
           ...
           'stepfunction2_mut10_IBM_1e4',...
           'stepfunction2_mut10_IBM_1e5',...
           'stepfunction2_mut10_IBM_1e6',...
           'stepfunction2_mut10_IBM_1e7',...
           ...
           'stepfunction3_mut10_IBM_1e4',...
           'stepfunction3_mut10_IBM_1e5',...
           'stepfunction3_mut10_IBM_1e6',...
           'stepfunction3_mut10_IBM_1e7',...
           ...
           'stepfunction4_mut10_IBM_1e4',...
           'stepfunction4_mut10_IBM_1e5',...
           'stepfunction4_mut10_IBM_1e6',...
           'stepfunction4_mut10_IBM_1e7',...
           ...
           'stepfunction5_mut10_IBM_1e4',...
           'stepfunction5_mut10_IBM_1e5',...
           'stepfunction5_mut10_IBM_1e6',...
           'stepfunction5_mut10_IBM_1e7',...
           };

% rundirs = {'constant_mut01_IBM_1e4',...
%            'constant_mut01_IBM_1e5',...
%            'constant_mut01_IBM_1e6',...
%            'constant_mut01_IBM_1e7',...
%            ...
%            'squarewave_mut01_IBM_1e4',...
%            'squarewave_mut01_IBM_1e5',...
%            'squarewave_mut01_IBM_1e6',...
%            'squarewave_mut01_IBM_1e7',...
%            ...
%            'sinusoidal_mut01_IBM_1e4',...
%            'sinusoidal_mut01_IBM_1e5',...
%            'sinusoidal_mut01_IBM_1e6',...
%            'sinusoidal_mut01_IBM_1e7',...
%            ...
%            'constant_mut10_IBM_1e4',...
%            'constant_mut10_IBM_1e5',...
%            'constant_mut10_IBM_1e6',...
%            'constant_mut10_IBM_1e7',...
%            ...
%            'squarewave_mut10_IBM_1e4',...
%            'squarewave_mut10_IBM_1e5',...
%            'squarewave_mut10_IBM_1e6',...
%            'squarewave_mut10_IBM_1e7',...
%            ...
%            'sinusoidal_mut10_IBM_1e4',...
%            'sinusoidal_mut10_IBM_1e5',...
%            'sinusoidal_mut10_IBM_1e6',...
%            'sinusoidal_mut10_IBM_1e7',...
%            };

clr=[187  59  14;
     221 118  49
     112 129  96
     216 197 147]./255;
 
 

outpath = '~/GitHub/Plankton_IBM/output/';

for i=1:numel(rundirs)
    
    outdir = [outpath rundirs{i}];
    
    cd(outdir)
    
    fntsz = 14;
    
    savefigs=true;
    
    
    IBM_output=load('IBM_output.mat');
    
    eco_params  = 
    .eco_params;
    env_forcing = IBM_output.env_forcing;
    genome      = IBM_output.bioinformatics.genome;
    rgb         = IBM_output.bioinformatics.rgb;
    
    t_opt       = eco_params.T_optimum;
    
    bio_out=IBM_output.yout(:,2:end-1);
    perc=(100.*bio_out./sum(bio_out,2));
    ycoord = linspace(eco_params.T_bins(1),eco_params.T_bins(end),size(perc,2));
    
    %% sample terminal leaf nodes
    nsample = 100; % max nsample terminal leaf nodes
    iextant = find(IBM_output.y(2:end-1));
    nextant = numel(iextant);
    nsample = min(nsample,nextant); % decrease if not enough nodes
    
    for j=1%:50
        isample = sort(randsample(iextant,nsample)); % sample without replacement
        
        %% bioinformatics toolbox
        
        
        genome_string = bioinf_fcns.print_genomes(genome(isample,:));
        
        
        % generate hamming-distance matrix
        L = eco_params.ngenome.*64;
        p = pdist(genome_string-'0','Hamming');
        p = squareform(p); % make square matrix
        X = p.*L./ eco_params.pneutral./2; % convert from fraction of bits different to total bits different
        % also account for mutation probability
        % 2-base jukes-cantor correction
        pcorr=p;
        pcorr(pcorr>0.5)=NaN;
        jk = -1/2 * log(1 - pcorr * 2/1);
        
        
        clear bingenes
        genmat=genome_string-'0';
        for k=1:nsample
            bingenes(k).Sequence=char(97+genmat(k,:)-'0');
            bingenes(k).Header='';
        end
        
        
        UPGMAtree = seqlinkage(jk.*L,'UPGMA',bingenes);
        NJtree    = seqneighjoin(jk.*L,'equivar',bingenes);
        
        UPGMAtree = reorder(UPGMAtree,1:nsample,'approximate',true);
        NJtree = reorder(NJtree,1:nsample,'approximate',true);
        
        h1 = plot(UPGMAtree,'orient','left');
        h2 = plot(NJtree,'orient','left');
        
        Xdist=[h1.LeafDots.XData h1.BranchDots.XData];
        Ydist=[h1.LeafDots.YData h1.BranchDots.YData];
        % Xdist=[h2.LeafDots.XData h2.BranchDots.XData];
        % Ydist=[h2.LeafDots.YData h2.BranchDots.YData];
        
        close(1)
        close(2)
        %%
        
        figure(111)
        subplot(2,3,ceil(i/4))
        hold on
        xx=sort(Xdist);
        xx=xx-max(xx);
        xx=xx(1:nsample);
        
        iclr = rem(i-1,4)+1;
        plot(xx,1:nsample,'color',[clr(iclr,:) 0.1],'LineWidth',2)
        set(gca,'XScale','lin')
        xlim([-1e3 -1])
    end
    
    xl=xlim;
    
    %% Coalescent theory
    
    for k=nsample:-1:2
        % DODGY IF MORE THAN 1 COALESCENCE PER GEN
        p  = 1./(2.*nextant.*(1./(k-1)-1./k));
        
        Tk_mn(nsample-k+1) = 1./p; % expected waiting time from coalescence k-1 to coalescence k
        Tk_sd(nsample-k+1) = sqrt((1-p)./p.^2);
    end
    
    plot(-[0 cumsum(Tk_mn)],nsample:-1:1,'-','color',clr(iclr,:).*0.75,'LineWidth',2)
    hold on
    % cumulative error calculated as root mean cumsum of squares
    plot(-[0 cumsum(Tk_mn)+sqrt(cumsum(Tk_sd.^2))],nsample:-1:1,'-','color',clr(iclr,:).*0.75,'LineWidth',1)
    plot(-[0 cumsum(Tk_mn)-sqrt(cumsum(Tk_sd.^2))],nsample:-1:1,'-','color',clr(iclr,:).*0.75,'LineWidth',1)
    % herrorbar(cumsum(Tk_mn),n-1:-1:1,sqrt(cumsum(Tk_sd.^2)),'m')
    
    box on
    
    drawnow
    
    
    
end
%%
ttls={'\DeltaT=0^\circC','\DeltaT=1^\circC','\DeltaT=2^\circC',...
      '\DeltaT=3^\circC','\DeltaT=4^\circC','\DeltaT=5^\circC'}
% ttls={'constant, \sigma_M = 0.01','diurnal, \sigma_M = 0.01','annual, \sigma_M = 0.01',...
%       'constant, \sigma_M = 0.10','diurnal, \sigma_M = 0.10','annual, \sigma_M = 0.10'}



for i=1:6
    h1=subplot(2,3,i);
    if i<=3
        h1.Position(2) = 0.525;
    else
        xlabel('Generations before end of simulation')
    end
    if i==1 |i==4
        ylabel('Lineages')
    end
    h1.Position(1) = 0.1+rem(i-1,3).*0.25;
    title(['(' char(96+i) ') ' ttls{i}])    
end

set(gcf,'color','w')
cd(outpath)
if savefigs
    export_fig -r450 coalesence.png
    print -depsc -tiff -r300 -painters coalesence.eps
end

%%





























return
%%
clear
clf

N=5000;
n=100;

for k=n:-1:2
    p  = 1./(2.*N.*(1./(k-1)-1./k));
    
    Tk_mn(n-k+1) = 1./p; % expected waiting time from coalescence k-1 to coalescence k
    Tk_sd(n-k+1) = sqrt((1-p)./p.^2);
end

plot(cumsum(Tk_mn),n-1:-1:1,'m-','LineWidth',2)
hold on
% cumulative error calculated as root mean cumsum of squares
plot(cumsum(Tk_mn)+sqrt(cumsum(Tk_sd.^2)),n-1:-1:1,'m-','LineWidth',1)
plot(cumsum(Tk_mn)-sqrt(cumsum(Tk_sd.^2)),n-1:-1:1,'m-','LineWidth',1)
% herrorbar(cumsum(Tk_mn),n-1:-1:1,sqrt(cumsum(Tk_sd.^2)),'m')


xlim([0 N])

N=5.*logspace(2,5,4)';
hold on
xx=linspace(0,max(Tk_mn),1000);
thalf = 2.*N./n;
plot(xx,n./(1+xx./thalf),'k--','LineWidth',2);




%%

% ET=0;
% for k=n:-1:2
%     Tk = 2.*N.*(1./(k-1)-1./k); % expected waiting time from coalescence k-1 to coalescence k
%
%     ET = ET+Tk;
%     plot(ET,k-1,'mo')
%     hold on
% end