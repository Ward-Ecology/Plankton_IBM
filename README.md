# Plankton_IBM

Instructions for reproducing the figures presented in "Rapid evolution allows coexistence of highly divergent lineages within the same niche", by Ward and Collins.

Published model results can be reproduced for different configurations by running the ``` Run_ecoevo.m ``` script. This will run the model and plot a selection of different figures. (Note that some simulations may take several hours or days to complete and may produce over 500 Gb of data.)

``` Run_ecoevo.m ``` includes the following default configuration:

``` Matlab
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
```

The environmental forcing function can be changed between the six different configurations shown in Figure 6 by changing the ``` eco_params.forcing ``` as follows:

| ``` eco_params.forcing ```  | Label in Figure 6 |
| ------------- | ------------- |
| ```constant```  | (a) constant temperature  |
| ```sinusoidal```  | (b) annual cycle  |
| ```stepfunction```  | (c) abrupt 5°C change  |
| ```squarewave```  | (d) diel cycle (speciation)  |
| ```gradual```  | (e) warming |
| ```gradualsinusoidal```  | (f) warming + annual cycle  |

The 1,000 year simulation shown in Figure 2 requires the following configuration:

``` Matlab
%% Overwrite default parameters here

eco_params.use_bioinf     = true; % if true, enable neutral genomes
eco_params.bioinf_verbose = false; % if true, save details of *all* cell divisions

env_forcing.tmax          = 1000.*360;        % (N.B. 360 days per year by default)
eco_params.forcing        = 'constant';   % define temperature forcing function
% forcing cases 'constant', 'stepfunction', 'sinusoidal', 'squarewave', 'gradual', 'gradualsinusoidal'

eco_params.ngenome        = 50; % number of neutral binary genes (53 bases each)
eco_params.pneutral       = 0.1;  % Set to 0.1 for 1000 year runs to avoid saturation
eco_params.resting_stages = false;

eco_params.V              = 1e-4; % Volume of culture (1e-6 is 1 ml) (All runs in Ward & Collins at 1e-4)
eco_params.nsuper         = 1;    % number of individuals per super-individual

eco_params.initialtraits = 'single'; % define initial trait distribution
```

These changes extend the model run time to 1000 years (```env_forcing.tmax = 1000.*360;```), supress verbose lineage tracking (```eco_params.bioinf_verbose = false;```), and decrease the probability of mutations in the binary genome (```eco_params.pneutral = 0.1;```). This last change is to prevent the number of neutral mutations from exceeding the size of the binary genome.

Figure 2 requires the model to be run using Beckmann et al.'s "trait diffusion" approximation of the individual based model. This can be done by setting ```eco_params.model_type = 'MCM';``` on line 10 of ``` Run_ecoevo.m ``` (with all other parameters set as standard).

Figure A.1 can be reproduced by running the ```clock_test.m``` script in the ```Sandbox``` directory.

Figures A.2 and A.3 can be reproduced by running ```Run_ecoevo.m ``` with the standard configuration modified such that ```eco_params.V = 1e-6;``` and ```eco_params.forcing = 'squarewave';```.

Note that ```Plankton_IBM``` requires the Matlab ```bioinformatics_toolbox```, ```distrib_computing_toolbox```, and ```statistics_toolbox```
