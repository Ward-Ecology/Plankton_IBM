%% Ecosystem equations
function [V eco_params bioinformatics phylogeny] = MUSE_IBM_eqs(t,V,env_forcing,eco_params,bioinformatics,bioinf_fcns,phylogeny)
 
% get environmental forcing (i.e. temperature)
    T  = env_forcing.Tfunc(t);
    
% get state variables from V_vect
    N  = V(1); 
    P  = V(2:eco_params.jpmax+1); 
    D  = V(eco_params.jpmax+2:end);
    
% index of live agents
    extant  = find(P); % index all non-zero populations
    
% evaluate limitation factors
    Nconc = N./eco_params.V;% convert to concentration for use in F_N
    F_N = Nutrient_function(Nconc,eco_params);
    F_T = Temperature_function(T,eco_params);
    
% remame parameters in shorter form
    mu0         = eco_params.m0_ln2_doubling;
    b_0         = eco_params.b_0;
    gamma       = eco_params.gamma_0_mort;
    tau0        = eco_params.tau_remin;
    nsuper      = eco_params.nsuper;
    T_optimum   = eco_params.T_optimum;
    delta_time  = eco_params.delta_t;
    gamma_crit  = gamma.*delta_time;
    
    nrgb = eco_params.nrgb;
    
    
% individual-based model equations
    %==========================
    
    % cell growth
    mu = mu0 .* b_0 .* F_T .* F_N .* (P>0);
    
    if eco_params.resting_stages
        isresting = find(F_T<eco_params.resting_cutoff);
        mu(isresting)    = 0.01 .* b_0 .* (P(isresting)>0);
        gamma_crit = gamma_crit .* ones(size(mu));
        gamma_crit(isresting) = 0.01.*delta_time;
        gamma_crit  = gamma_crit(extant);
    end
    
    % stochastic mortality
    i_dying = extant(find(rand(nnz(extant),1) < gamma_crit));
    Pdead   = P(i_dying); % sum of dying cells to be passed to detritus
    % apply deaths
    P(i_dying)         = 0; % kill cells
    T_optimum(i_dying) = -999; % and set T_optimum to -999
    mu(i_dying) = 0; % cells that die this timestep do not grow
    
    % index cells with sufficient biomass to divide
    i_dividing = find(P > 2.*b_0);
    
    % rates of change for continuous state variables (N and D)        
    dNdt = - sum(mu .* nsuper) + tau0.*D; % N uptake = sum of cell growth * number of cells per agent
    dDdt = - tau0.*D; % detrital remineralisation
    
    % halve biomass of dividing cells
    P(i_dividing) = P(i_dividing)./2;
    extinct = find(~P);          % identify vacant spaces in P vector
    numdiv  = numel(i_dividing); % check number of vacant slots for dividing c
    
    if numdiv>numel(extinct) % if there is not sufficient space for cells
        % error if P vector is too small
        error('agent overflow (agent vector not big enough)');
    end
    
    tofill     = extinct(1:numdiv); % get index for new cells
    P(tofill)  = P(i_dividing);     % duplicate dividing cells
    % mutate dividing cells
    T_optimum(tofill)     = T_optimum(i_dividing) ...
                          + randn(size(i_dividing)).*eco_params.sigmaM; % duplicate thermal optima (mutations at start of timestep)
    
    if eco_params.use_bioinf
        % get time cells die or divide
        if eco_params.bioinf_verbose
            ind_die = sort([phylogeny.phylogeny_index(i_dying); phylogeny.phylogeny_index(i_dividing)]);
            phylogeny.time_of_dying(eco_params.nwritten(3)+(1:size(ind_die,1)),:) = [ind_die ones(size(ind_die)).*t];
            eco_params.nwritten(3) = eco_params.nwritten(3) + size(ind_die,1);
            phylogeny.phylogeny_index(i_dying)      = 0; % delete dying cells from phylogeny_index
        end
        
        % update extant agent indices
        bioinformatics.rgb(i_dying,:)           = 0;
        bioinformatics.otu(i_dying,1)           = 0;
        bioinformatics.genome(i_dying,:)        = 0;
        
        n_i_div = numel(i_dividing); % number of dividing agents
        if n_i_div>0
            
            % copy genomes and taxonomic ID
            bioinformatics.rgb(tofill,:)    = bioinformatics.rgb(i_dividing,:);
            bioinformatics.genome(tofill,:) = bioinformatics.genome(i_dividing,:);
            bioinformatics.otu(tofill,1)    = bioinformatics.otu(i_dividing,1);
            
            % rgb mutates every timestep
            i_mutates = [tofill i_dividing];
            bioinformatics.rgb(i_mutates,:)    = bioinf_fcns.rgb_mutate(    bioinformatics.rgb(i_mutates,:),ones(size(i_mutates)));
            % binary genes may have lower probability
            mutrand = rand(n_i_div,1);
            imut    = find(mutrand<=eco_params.pneutral);
            if ~isempty(imut)
                i_mutates = [tofill(imut) i_dividing(imut)];
                bioinformatics.genome(i_mutates,:) = bioinf_fcns.gene_mutate(bioinformatics.genome(i_mutates,:),ones(size(i_mutates)), eco_params.pmut);
            end
            
            if eco_params.bioinf_verbose
                % get index of dividing cell
                ind_old = phylogeny.phylogeny_index(i_dividing);
                % get index of new daughter cell
                ind_new = eco_params.nwritten(1)+(1:2:2.*numel(i_dividing))';
                
                % keep track of written divisions
                eco_params.nwritten(1) = max(ind_new)+1;
                n_offset               = eco_params.nwritten(2);
            
                iii=i_dividing;
                for i=0:1
                    % write parent ID for new cell
                    phylogeny.phylodata(ind_new-n_offset+i,1) = ind_old;
                    % 
                    phylogeny.phylodata(ind_new-n_offset+i,2) = ind_new+i;
                    phylogeny.phylodata(ind_new-n_offset+i,3) = t;
                    phylogeny.phylodata(ind_new-n_offset+i,4) = T_optimum(iii);
                    phylogeny.phylodata(ind_new-n_offset+i,5) = bioinformatics.otu(iii,:);
                    phylogeny.phylodata(ind_new-n_offset+i,6:(5+eco_params.nrgb))   = bioinformatics.rgb(iii,:);
                    if eco_params.bioinf_verbose
                        phylogeny.phylodata(ind_new-n_offset+i,(6+eco_params.nrgb):end) = bioinformatics.genome(iii,:);
                    end
                    
                    % update indices
                    phylogeny.phylogeny_index(iii)  = ind_new+i;   % change node index of dividing cell
                    iii = tofill;
                end
            end
            
        end
        if nnz(nnz(P))==0
            error('No live agents remaining.')
        end
    end
    
    % ecosystem dynamics 
    N = N + dNdt .* delta_time; 
    P = P + mu   .* delta_time;
    D = D + dDdt .* delta_time + sum(Pdead .* nsuper);
    
    % Update traits
    eco_params.T_optimum = T_optimum;
   
% Collate ouptut vector
    V=[N;P;D];

end