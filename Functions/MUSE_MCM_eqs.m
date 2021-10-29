%% Ecosystem equations
function [DV,mu] = MUSE_MCM_eqs(t,V0,env_forcing,eco_params)

% get environmental forcing (i.e. temperature)
    T  = env_forcing.Tfunc(t);
    
    
% get state variables from V_vect
    N  = V0(1);
    P  = V0(2:eco_params.jpmax+1); 
    D  = V0(eco_params.jpmax+2:end);
    
% evaluate limitation factors
    F_T = Temperature_function(T,eco_params)';
    F_N = Nutrient_function(N,eco_params);
    
% remame parameters in shorter form
    mu0    = eco_params.mu0_max_rate;
    gamma  = eco_params.gamma_0_mort .* ones(size(P));
    tau0   = eco_params.tau_remin;
    
% ecosystem model equations
    %==========================
    % ODEs
    
    mu  = mu0.*F_T.*F_N;
    
    if eco_params.resting_stages
        isresting = find(F_T<eco_params.resting_cutoff);
        mu(isresting)    = 0.01;
        gamma(isresting) = 0.01;
    end
    
    muP  = mu.*P;
    Tmut = eco_params.TT.*mu';
    
    dPdt = + muP + Tmut*P - gamma .* P;
    dNdt = - sum(muP) + tau0 .* D;
    dDdt = + sum(gamma .* P) - tau0 .* D;
   
    
% Collate ouptut vector
    DV=[dNdt; dPdt; dDdt];

end