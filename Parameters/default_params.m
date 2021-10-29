% load function libraries
forcing_fcns = forcing_functions; 
bioinf_fcns  = bioinf_functions; 

%% Set default parameters
%% General parameters
eco_params.R_mass_conc     = 5;                 % input nutrient concentration     (mmol N m^-3)
eco_params.kN_half_sat     = 0.15;              % nutrient half-saturation         (mmol N m^-3)
eco_params.gamma_0_mort    = 0.1;               % mortality rate                   (d^-1)
eco_params.tau_remin       = 0.25;              % remineralization rate            (d^-1)
eco_params.theta_norm      = 6;                 % width of thermal reaction norm   (degrees C)
eco_params.mu0             = log(2);            % maximum specific growth rate

%% Individual-based model
eco_params.V               = 1e-4;              % Volume of culture                (1 ml)
eco_params.b_0             = 5e-10;             % Minimum individual biomass       (mmol N)
eco_params.sigmaM          = 0.1;               % standard deviation of mutations  (degrees C)
eco_params.nsuper          = 1;                 % # cells per super individual
eco_params.n_agents_init   = 1e5;               % initial # agents 
eco_params.delta_t         = 1/24;              % time step                        (days)

%% Multi-compartment model
eco_params.deltaT_opt      = 0.1;               % temperature width of genotypes   (degrees C)
eco_params.delta           = 1/3;               % genotype transfer factor         (-)

%% Sensitivity experiments
eco_params.resting_stages  = false;             % resting stages not set by default
eco_params.resting_cutoff  = 0.1;
%% Output
eco_params.nsample         = 1e2;               % samples of IBM at each output timestep

%% Trait space
eco_params.initialtraits   = 'single';          % one initial plankton cell
eco_params.mintemp         =  5;                % minimum thermal optimum in MCM
eco_params.maxtemp         =  25;               % maximum thermal optimum in MCM

%% set environmental forcing (set in forcing_functions)
eco_params.forcing        = 'sinusoidal';       
env_forcing.daysperyear   = 360;  % (N.B. 360 days per year)
env_forcing.tmax          = 10*env_forcing.daysperyear;  % (N.B. 360 days per year)
env_forcing.period        = 1*env_forcing.daysperyear;   % (days)

%% bioinformatics
eco_params.use_bioinf       = true;
eco_params.bioinf_verbose   = false; % apply cautiously!
                                     % can generate VERY large output files
                                     % i.e. 10s to 100s of Gb
eco_params.nrgb             = 3;    % number of rgb genes
eco_params.ngenome          = 0;    % number of binary genes
eco_params.pneutral         = 1;    % probability of a binary gene mutation