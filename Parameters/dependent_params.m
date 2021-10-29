%% set environmental temperature as function of time
switch eco_params.forcing
    case 'constant'
        env_forcing.Tfunc = forcing_fcns.constant;
    case 'stepfunction'
        env_forcing.Tfunc = forcing_fcns.stepfunction;
    case 'sinusoidal'
        env_forcing.Tfunc = forcing_fcns.sinusoidal;
    case 'gibbs_effect'
        env_forcing.Tfunc = forcing_fcns.gibbs_effect;
    case 'sinustep'
        env_forcing.Tfunc = forcing_fcns.sinustep;
    case 'squarewave'
        env_forcing.Tfunc = forcing_fcns.squarewave;
    case 'gradual'
        env_forcing.Tfunc = forcing_fcns.gradual;
    case 'gradualsinusoidal'
        env_forcing.Tfunc = forcing_fcns.gradualsinusoidal;
end


%%
switch eco_params.model_type
%%
    case 'IBM'
        disp('Initialising Individual-Based Model') 
        
        % set maximum doubling rate from eco_params.mu0
        eco_params.m0_ln2_doubling = eco_params.mu0./log(2); % maximum doubling rate (d^-1)
        
        % calculate the maximum possible number of agents
        % (total biomass conc * culture volume / minimum individual biomass / individuals per agent)
        eco_params.jpmax = round(eco_params.R_mass_conc .* eco_params.V ...
                        ./ eco_params.b_0 ./ eco_params.nsuper); 
        
        % vector of thermal optima (uniform initial distribution)
        eco_params.T_optimum = zeros(eco_params.jpmax,1);
        
        if strmatch(eco_params.initialtraits,'uniform') 
            % uniform initial trait distribution
            rng(1); % set repeatable random seed
            eco_params.T_optimum(1:eco_params.n_agents_init) = eco_params.mintemp ...
                + (eco_params.maxtemp-eco_params.mintemp) .* sort(rand(eco_params.n_agents_init,1));
        elseif strmatch(eco_params.initialtraits,'gaussian')
            % gaussian initial trait  distribution (equivalent to Beckmann's constant-forcing steady-state)
            rng(1); % set repeatable random seed
            eco_params.T_optimum(1:eco_params.n_agents_init) = 15 ...
                                                             + sort(randn(eco_params.n_agents_init,1).*0.855);
        elseif strmatch(eco_params.initialtraits,'single')
            % Initial population is a single cell
            eco_params.n_agents_init = 1;
            eco_params.T_optimum(1)  = 15;
        end
        % reset random seed
        rng('shuffle')
        
        % Define thermal optima bins for output visualisation
        eco_params.T_bins    = [eco_params.mintemp-eco_params.deltaT_opt/2:...
                                eco_params.deltaT_opt:...
                                eco_params.maxtemp+eco_params.deltaT_opt/2];
        
        % initialise agent vector to maximum possible size
        P0 = zeros(eco_params.jpmax,1);
        
        % The phytoplankton population is initialized by randomly assigning
        % the biomass of each initial agent (distributed uniformly between minimum 
        % bo and maximum 2*bo biomass)
        P0(1:eco_params.n_agents_init) = eco_params.b_0 ...
                                       + eco_params.b_0 .* rand(eco_params.n_agents_init,1);
        
        % assign remaining N budget (R) to N
        N0    = eco_params.R_mass_conc .* eco_params.V - sum(P0).*eco_params.nsuper;
        if N0<0
            error('Initial population too high!\n%s','(Sum(P0)>eco_params.R_mass_conc)')
        end
        
        % set initial detritus to zero
        D0    = 0;
        
        % initialise state variables as vector
        y0 = [N0;P0;D0];
        
        % preallocate
        Topt_out = zeros(env_forcing.tmax,eco_params.nsample).*NaN;
%%
    case 'MCM'
        disp(' Initialising Multi-Compartment Model')
        
        % set maximum growth rate from eco_params.mu0
        eco_params.mu0_max_rate = eco_params.mu0; % maximum growth rate (d^-1)
        
        % vector of thermal optima
        eco_params.T_optimum = eco_params.mintemp:eco_params.deltaT_opt:eco_params.maxtemp;
        
        % number of plankton populations
        eco_params.jpmax     = numel(eco_params.T_optimum);
        eco_params.T_bins    = [eco_params.mintemp:eco_params.deltaT_opt:eco_params.maxtemp];
        eco_params.delta     = eco_params.sigmaM^2 / (3*eco_params.deltaT_opt^2);

        eco_params.deltaM = eco_params.sigmaM.^2 ...      % genotype transfer function
            ./ ( 3.*eco_params.deltaT_opt.^2 );
        n      = eco_params.jpmax;                      % n 'genotypes'
        delta  = ones(n,1).*eco_params.delta;          % expanded vector of genotype transfer function
        TT     = spdiags([delta delta],[-1 1], n, n); % set off-diagonals to delta M
        TT     = spdiags(-sum(TT,2),0, TT);             % set diagonal to conserve mass
%         TT = TT + speye(size(TT));
        
        eco_params.TT=TT;
        
        % set initial conditions
        % set blank vector with size [1,jpmax]
        P0    = zeros(eco_params.jpmax,1);
        % set non-zero initial biomass for some population/s
        init_bio = eco_params.b_0 + eco_params.b_0 .* rand(eco_params.n_agents_init,1);
%         P0(:) = sum(init_bio).*eco_params.nsuper./eco_params.jpmax; % initial uniform across T_opt biomass
        
        P0(round(size(P0,1)/2)) = sum(init_bio).*eco_params.nsuper./eco_params.jpmax;
        
        % assign remaining N budget (R) to N
        N0    = eco_params.R_mass_conc - sum(P0);
        % set initial detritus to zero
        D0    = 0;
        
        % initialise state variables as vector
        y0 = [N0;P0;D0];
        
        % preallocate
        bio_out = zeros(env_forcing.tmax,eco_params.jpmax);
        eco_params.nsample = eco_params.jpmax;
    otherwise
        %%
        disp('Unknown model type (select ''IBM'' or ''MCM'').')
end



% Bioinformatics
if eco_params.use_bioinf
    % each population (MCM) or individual (IBM) can be assigned an rgb 'genome' of eco_params.nrgb genes
    % each 'genome' is a vector of signed (i.e. +/-) floating point numbers initialised at [0,0,...,0]
    % the genome is ecologically neutral, with no effect on ecological dynamics
    % mutations are equivalent to a random walk in an n-dimensional space
    % (where n = eco_params.nrgb)
    % mutations are carried out by bioinf_fcns.rgb_mutate
    rgb        = zeros(eco_params.jpmax,eco_params.nrgb);
    
    % fixed taxonomic ID assigned at start of run
    otu        = zeros(eco_params.jpmax,1);
    otu(1:eco_params.n_agents_init,1)=1:eco_params.n_agents_init;

    % populations/individuals can also be assigned a binary genome of 64*eco_params.ngenome bits
    % the binary genome is also ecologically neutral, with no effect on ecological dynamics
    % mutations occur as random bit flips, with one randomly selected bit
    % selected to flip at each generation for each population/individual
    % mutations are carried out by bioinf_fcns.gene_mutate
    %
    %
    % to decrease memory load, the 53 bit genes are stored as decimal integers
    % this way, 53 ones and zeros (a character string of 106 bytes) 
    % can be stored as a single decimal digit (double precision scalar of 8 bytes)
    % 000...000 --> 2^0 -1 --> 0
    % 111...111 --> 2^53-1 --> 9.0072e+15
    %
    % for spase arrays
    % N.B. 53 is the maximum length that can be safely converted between binary and decimal format
    % bin2dec(dec2bin(2^53-1)) == 2^53-1 --> true
    % bin2dec(dec2bin(2^53)) --> Error (Binary character vector must be 53 bits or less.)
    % for full arrays, use uint64 - allowing 64 bits per gene
    genome     = sparse(eco_params.jpmax,eco_params.ngenome);
    genome     = uint64(full(genome));

    % Set up structural array for bioinformatoc data
    phylogeny=struct([]); 
    switch eco_params.model_type
        case 'IBM'
            % Clear output directory
            !rm ../output/*
            if eco_params.bioinf_verbose
                clear phylogeny
                % seed individuals assigned srcnode ID of 1
                srcnode = ones(eco_params.n_agents_init,1);
                % seed individuals assigned consecutive snknode IDs of 1:eco_params.n_agents_init
                snknode = (1:eco_params.n_agents_init)';
                % seed individuals assigned division (i.e. birth) time of zero
                divtime = zeros(eco_params.n_agents_init,1);
                % thermal optima recorded
                t_optim = eco_params.T_optimum(1:eco_params.n_agents_init);
                % seed individuals assigned rgb gene of [0 0 0]
                rgb_phy = zeros(eco_params.n_agents_init,eco_params.nrgb);
                % seed individuals assigned consecutive heritable OTU numbers
                otu_phy = otu(1:eco_params.n_agents_init);
                if eco_params.bioinf_verbose
                    % generate array to hold binary genomes
                    gen_phy = zeros(eco_params.n_agents_init,eco_params.ngenome);
                else
                    gen_phy = zeros(eco_params.n_agents_init,0);
                end
                % assign index numbers for seed individuals
                phylogeny_index = zeros(eco_params.jpmax,1);
                phylogeny_index(1:eco_params.n_agents_init)=1:eco_params.n_agents_init;

                % initialise output .csv files                
                phylodata = [srcnode,snknode,divtime,t_optim,otu_phy,rgb_phy,gen_phy];
                eco_params.fmt = ['%i,%i,%.4f,%.4f,%i' repmat(',%.4f',1,eco_params.nrgb) repmat(',%u',1,size(gen_phy,2)) '\n'];
                [fid,msg] = fopen('../output/phylogeny_0001.csv', 'w');
                fprintf(fid,eco_params.fmt,phylodata');
                fclose(fid);

                % note n of individuals written (keeps track of output file size)
                eco_params.nwritten = [eco_params.n_agents_init eco_params.n_agents_init];

                % reset to all zeros
                phylogeny.phylodata       = zeros(5e6,size([srcnode,snknode,divtime,t_optim,otu_phy,rgb_phy,gen_phy],2));
                phylogeny.phylogeny_index = phylogeny_index;
                phylogeny.time_of_dying   = zeros(1e7,2);
                eco_params.nwritten(3)    = 0;
            end
    end
    
    % create structiral arrays for concise passing to functions
    bioinformatics.rgb    = rgb;
    bioinformatics.otu    = otu;
    bioinformatics.genome = genome;
    
    clear rgb genome phylodata phylogeny_index 
else 
    bioinformatics=[];
    phylogeny=[];
end

%% make mutation matrix




