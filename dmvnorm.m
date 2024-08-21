% SW: This method is private. It is the same as mvtnorm::dmvnorm() function, 
% but without sanity checks for sigma. We perform the sanity checks before.
function logretval = dmvnorm(x, mean, sigma, loga)
% calculates the multivariate normal pdf 
% input
    % x (m x n) with m the number of query points and n the dimensionality of mean
    % mean (nx1) 
    % sigma (nxn)
    % loga optional flag to take the logarithmic of the pdf
    
    if nargin < 4
        loga = false;
    end

    logretval = mvnpdf(x, mean', sigma );
    if loga
        logretval = log(logretval);
    end
end