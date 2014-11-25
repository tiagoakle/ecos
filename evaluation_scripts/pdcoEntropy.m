function [f,g,H] = pdcoEntropy(x)
    f = sum(x.*log(x));
    g = log(x)+1;
    H = diag(sparse(1./x));
end

