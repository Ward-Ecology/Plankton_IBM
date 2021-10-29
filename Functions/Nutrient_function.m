%% Nutrient_function
function [F_N] = Nutrient_function(N,eco_params)
    
    % remame parameters in shorter form
    kN = eco_params.kN_half_sat;
    
    % caculate function
    F_N = N ./ (kN + N); % Beckmann equation 9
end