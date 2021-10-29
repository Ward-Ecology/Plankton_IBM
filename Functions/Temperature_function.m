%% Temperature_function
function [F_T] = Temperature_function(T,eco_params)
    
    % remame parameters in shorter form
    T_opt = eco_params.T_optimum;
    theta = eco_params.theta_norm;
    
    % calculate function
    F_T = exp( -((T - T_opt)./theta).^2 ); % Beckmann equation 8
    % Note that this function is a significant (and unavoidable?) computational bottleneck
    
end