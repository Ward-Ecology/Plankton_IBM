function [ forcing_fcns ] = forcing_functions ( )
    forcing_fcns.constant           = @constant;
    forcing_fcns.stepfunction       = @stepfunction;
    forcing_fcns.sinusoidal         = @sinusoidal;
    forcing_fcns.gibbs_effect       = @gibbs_effect;
    forcing_fcns.gradual            = @gradual;
    forcing_fcns.sinustep           = @sinustep;
    forcing_fcns.squarewave         = @squarewave;
    forcing_fcns.gradualsinusoidal  = @gradualsinusoidal;
end

%%
function T = constant(t)
    T = ones(size(t)).*15;
end

%%
function T = stepfunction(t)
    T(t> 50.*360) = 15+3;
    T(t<=50.*360) = 15;
end

%%
function T = sinusoidal(t)
    T = 15 + 5 * sind(t-60); % -60 offsets sine wave so population
                              % is more evenly balanced between hot and cold at end of years
       
%     T(t<1*360-120) =15;
end


%%
function T = gibbs_effect(t)                             
               
    t60=t-60;   % -60 offsets sine wave so population
                % is more evenly balanced between hot and cold at end of years
    T = 15 + 5 * (sind(t60)+sind(3.*t60)./3+sind(5.*t60)./5);

    T(t<1*360-120) =15;
end

%%
function T = sinustep(t)
    T(t> 10.*360) = 20;
    T(t<=10.*360) = 15;

    T = T + 5 * sind(t); % -60 offsets sine wave so population
                              % is more evenly balanced between hot and cold at end of years
end


%%
function T = squarewave(t)

    T = zeros(size(t));
    
    t_dec=rem(t,1);
    T(t_dec< 0.5)   = 10;
    T(t_dec>=0.5)   = 20;
    
    T(t   <=5.*360) = 15; 

end



%%
function T = gradual(t)

    T = 15 + 1/3.*(t-5*360)./360;
    
    T(t<=5*360) =15;
    
end

%%
function T = gradualsinusoidal(t)

    yr= t./360;
    
    T = 15 + (yr-5)./2; % (1 degree per 20 years)
    T(yr<5) =15;    
    
    T = T + 5 * sind(t); 
end




