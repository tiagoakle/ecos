%Returns a handle to a function that calculates 
% the value gradient and objective for 
% sum_i x_ic_i^{-1}\log(x_ic_i^{-1})

function func_handle =  get_pdcoEntropy(c);
   cinv = 1./c;
   func_handle = @(x)(pdcoEntropyScaled(x,cinv));
end

function [f,g,H] = pdcoEntropyScaled(x,s)
    xs = x.*s;
    f = sum(xs.*log(xs));
    g = s.*(log(xs)+1);
    H = diag(s.*sparse(1./x));
end

