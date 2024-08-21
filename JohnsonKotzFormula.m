function result = JohnsonKotzFormula(mean, sigma, lower, upper)
% Compute truncated mean and truncated variance when only a subset of 
% k < n x_1,..,x_k variables are truncated, using Johnson/Kotz (1972, p.70)
% formula.
% input
    % mean Mean vector (k x 1)
    % sigma Covariance matrix (k x k)
    % lower Lower truncation point (k x 1)
    % upper Upper truncation point (k x 1)
% output
    % structure with fields tmean and tvar containing mean and variance of
    % truncated distribution.

    % Determine which variables are truncated
    idx = find(~isinf(lower) | ~isinf(upper)); % Index of truncated variables
    n = length(mean);
    k = length(idx); % Number of truncated variables
    if k >= n
        error(sprintf("Number of truncated variables (%d) must be lower than total number of variables (%d).", k, n));
    end
    if k == 0
        result = struct('tmean', mean, 'tvar', sigma); % No truncation
        return;
    end

    % Transform to zero mean first
    lower = lower - mean;
    upper = upper - mean;

    % Partitioning of sigma
    % sigma = [ V11  V12 ]
    %         [ V21  V22 ]
    V11 = sigma(idx, idx);
    V12 = sigma(idx, setdiff(1:n, idx));
    V21 = sigma(setdiff(1:n, idx), idx);
    V22 = sigma(setdiff(1:n, idx), setdiff(1:n, idx));

    % Determine truncated mean xi and truncated variance U11
    r = mtmvnorm(zeros(k, 1), V11, lower(idx), upper(idx));
    xi = r.tmean;
    U11 = r.tvar;

    invV11 = inv(V11); % V11^(-1)

    % See Johnson/Kotz (1972), p.70 formula
    tmean = zeros(n, 1);
    tmean(idx) = xi;
    tmean(setdiff(1:n, idx)) = xi' * invV11 * V12;
    tvar = nan(n, n);
    tvar(idx, idx) = U11;
    tvar(idx, setdiff(1:n, idx)) = U11 * invV11 * V12;
    tvar(setdiff(1:n, idx), idx) = V21 * invV11 * U11;
    tvar(setdiff(1:n, idx), setdiff(1:n, idx)) = V22 - V21 * (invV11 - invV11 * U11 * invV11) * V12;

    tmean = tmean + mean;

    result = struct('tmean', tmean, 'tvar', tvar);
end