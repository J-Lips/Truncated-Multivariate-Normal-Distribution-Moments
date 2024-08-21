function cargs = checkTmvArgs(mean, sigma, lower, upper)
% checks sanity of input arguments
% input
    % mean Mean vector (k x 1)
    % sigma Covariance matrix (k x k)
    % lower Lower truncation point (k x 1)
    % upper Upper truncation point (k x 1)
% output
    % structure with fields mean, sigma, lower, upper 

    if isempty(lower) || any(isnan(lower))
        error('lower not specified or contains NA');
    end
    if isempty(upper) || any(isnan(upper))
        error('upper not specified or contains NA');
    end
    if ~isnumeric(mean) || ~isvector(mean)
        error('mean is not a numeric vector');
    end
    if isempty(sigma) || any(any(isnan(sigma)))
        error('sigma not specified or contains NA');
    end
    
    if ~ismatrix(sigma)
        sigma = asmatrix(sigma);
    end
    
    if size(lower, 1) ~= size(upper, 1)
        error('lower and upper have non-conforming size');
    end
    
    checkSymmetricPositiveDefinite(sigma);
    
    d = length(mean);
    
    if length(mean) ~= size(sigma, 1)
        error('mean and sigma have non-conforming size');
    end
    
    if length(lower) ~= size(lower, 1) || length(upper) ~= size(upper, 1)
        error('lower and upper have non-conforming size');
    end
    
    cargs = struct('mean', mean, 'sigma', sigma, 'lower', lower, 'upper', upper);
end

function checkSymmetricPositiveDefinite(sigma)
% checks sanity of sigma 
% input
    % sigma Covariance matrix (k x k)
% output
    % -
    if ~isequal(sigma, sigma')
        error('sigma must be symmetric');
    end
    if any(eig(sigma) <= 0)
        error('sigma must be positive definite');
    end
end