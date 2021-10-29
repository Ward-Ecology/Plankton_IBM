function y = exp_lut(x);
% approximate output of (slow) exponential function from 
% precalculated look-up table with linear approximation.
    persistent xLUT yLUT
    if ( isempty(yLUT))
        xLUT = linspace(-20,45,1000);
        yLUT = exp(xLUT);
    end
    
    y = interp1(xLUT,yLUT,x);
    % set y for x under range to zero
    y(x<xLUT(1)) = 0;
end