% Plankton IBM based on Beckmann et al. 2019 - "Phytoplankton adaptation in ecosystem models"

clear,clc

addpath Plotting
addpath Functions
addpath Post_processing

% set model (IBM or MCM)
eco_params.model_type = 'IBM';

%% set ecophyslogical parameters
run('parameters/default_params.m');

%% Overwrite default parameters here

eco_params.use_bioinf     = true; % if true, enable neutral genomes
eco_params.bioinf_verbose = true; % if true, save details of *all* cell divisions

env_forcing.tmax          = 15.*360;        % (N.B. 360 days per year by default)
eco_params.forcing        = 'constant';   % define temperature forcing function
% forcing cases 'constant', 'stepfunction', 'sinusoidal', 'squarewave', 'gradual', 'gradualsinusoidal'

eco_params.ngenome        = 50; % number of neutral binary genes (53 bases each)
eco_params.pneutral       = 1;  % Set to 0.1 for 1000 year runs to avoid saturation
eco_params.resting_stages = false;

eco_params.V              = 1e-4; % Volume of culture (1e-6 is 1 ml) (All runs in Ward & Collins at 1e-4)
eco_params.nsuper         = 1;    % number of individuals per super-individual

eco_params.initialtraits = 'single'; % define initial trait distribution

%% calculate dependent parameters
run('parameters/dependent_params.m');

%% run simulation
y = y0;

switch eco_params.model_type
   
    case {'IBM'}  
        disp('   Evaluating Individual-Based Model') 
        [yout y tout eco_params bioinformatics phylogeny] = timestep_IBM(y,env_forcing,eco_params,bioinformatics,bioinf_fcns,phylogeny);
        % yout         = binned state variable output (all timesteps)
        % y            = individual state variable output (last timestep)
        % eco_params   = ecological parameters
        
%         profile viewer

        
        save('output/IBM_output.mat', '-v7.3')
        
        location = [pwd '/output/'];
        if eco_params.bioinf_verbose
            post_process_phylogeny(location)
            final_year_phylogeny(location)
            plot_IBM_figures
        else
            plot_IBM_nonverbose
        end

                        
    case{'MCM'}
        
        disp('Evaluating Multi-Compartment Model') 
        if true
            % use ode45 to solve MCM model
            % ODE45
            [tout,yout] = ode45(@(t,y)  MUSE_MCM_eqs(t,y,env_forcing,eco_params), ... % model function
                                        eco_params.delta_t:env_forcing.tmax,... % vector defining all data output times during defined integration period
                                        y0); % initial condition%% ODE45 - VariabT 1520C
        else
            % use forward Euler to solve MCM model
            % (N.B. IBM olnly runs with Forward Euler, so need to check this works as well as ode45 for MCM case)
            % Forward Euler
            yout   = zeros(env_forcing.tmax,size(T_bins,2)+2);
            for t=0:eco_params.delta_t:env_forcing.tmax
                [DV,~] = MUSE_MCM_eqs(t,y,env_forcing,eco_params);
                y=y+DV.*eco_params.delta_t;
                if rem(t,1)==0 & t>0
                    yout(t,:) = y;
                end
            end
            tout=1:env_forcing.tmax;
        end
        
        %% Offline diagnostics
        % ode45 cannot output diagnostic variables
        % here we re-evaluate MCM at set timepoints to get bioinformatic output
        disp('Re-evaluating Multi-Compartment Model (diagnostics)') 
        clear bio_out rgb_out genome_out
        bio_out    = zeros(env_forcing.tmax,eco_params.jpmax);
        rgb_out    = zeros(env_forcing.tmax,eco_params.jpmax,eco_params.nrgb);
        genome_out = zeros(env_forcing.tmax,eco_params.jpmax,eco_params.ngenome);
        genome_out = uint64(genome_out);
        
        rgb    = zeros(eco_params.jpmax,eco_params.nrgb);
        genome = uint64(zeros(eco_params.jpmax,eco_params.ngenome));
        
        mu_out = zeros(eco_params.jpmax,env_forcing.tmax);
        
        for dy=1:env_forcing.tmax
            biomass = yout(dy,2:end-1)';
            bio_out(dy,:)   = biomass;
            
            % call MUSE_MCM_eqs to get growth rate at t = tout(dy)
            [~,mu] = MUSE_MCM_eqs(tout(dy),yout(dy,:)',env_forcing,eco_params);
            
            mu_out(:,dy) = mu;
            
%             plot(mu)
%             axis([0 301 0 1.5])
%             drawnow
            
            if eco_params.use_bioinf
                    
                mutrand = rand(eco_params.jpmax,1); % vector of uniform random numbers: U[0,1]
                imut    = find(mu.*eco_params.pneutral>mutrand); % if specific growth rate > U[0,1]
                rgb(imut,:)    = bioinf_fcns.rgb_mutate(    rgb(imut,:),ones(size(imut)));
                genome(imut,:) = bioinf_fcns.gene_mutate(genome(imut,:),ones(size(imut)));

                % use muP to work out neutral genetic descent
                [jj] = bioinf_fcns.bioinf_inherit(biomass,mu,eco_params.TT);
                rgb     = rgb(jj,:);
                genome  = genome(jj,:);

                rgb_out(dy,:,:) = rgb;
                genome_ndays(dy,:,:) = genome;
            end
        end
        
        disp('Saving Multi-Compartment Model') 
        save('output/MCM_output.mat')
        
        disp('Plotting Multi-Compartment Model') 
        plot_MCM_figures
end





















