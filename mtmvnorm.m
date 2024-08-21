function result = mtmvnorm(mean, sigma, lower, upper, doComputeVariance)
% Mean and Covariance of the truncated multivariate distribution (double truncation, general sigma, general mean)
% input
    % mean Mean vector (k x 1)
    % sigma Covariance matrix (k x k)
    % lower Lower truncation point (k x 1)
    % upper Upper truncation point (k x 1)
    % doComputeVariance (optional) flag whether to compute variance
% output
    % structure with fields tmean and tvar containing mean and variance of
    % truncated distribution.
    
    if nargin < 5
        doComputeVariance = true;
    end
    N = length(mean);

    % Check input parameters
    cargs = checkTmvArgs(mean, sigma, lower, upper);
    mean = cargs.mean;
    sigma = cargs.sigma;
    lower = cargs.lower;
    upper = cargs.upper;

    % Check number of truncated variables; if only a subset of variables is truncated
    % we can use the Johnson/Kotz formula together with mtmvnorm

    % Determine which variables are truncated
    idx = find(~isinf(lower) | ~isinf(upper)); % Index of truncated variables
    k = length(idx); % Number of truncated variables
    if k < N
        result = JohnsonKotzFormula(mean, sigma, lower, upper);
        return;
    end

    % Truncated Mean
    TMEAN = zeros(N, 1);
    % Truncated Covariance matrix
    TVAR = nan(N, N);

    % Shift the integration limits by -mean to make the mean 0
    a = lower - mean;
    b = upper - mean;
    lower = lower - mean;
    upper = upper - mean;

    % One-dimensional marginal density
    F_a = zeros(N, 1);
    F_b = zeros(N, 1);

    zero_mean = zeros(N, 1);

    % Pre-calculate one-dimensional marginals F_a[q] once
    for q = 1:N
        tmp = dtmvnorm_marginal([a(q), b(q)], q, zero_mean, sigma, lower, upper, false);
        F_a(q) = tmp(1);
        F_b(q) = tmp(2);
    end

    % 1. Determine E[X_i] = mean + Sigma %*% (F_a - F_b)
    TMEAN = sigma * (F_a - F_b);

    if doComputeVariance
        % TODO:
        % Calculating the bivariate densities is not necessary
        % in case of conditional independence.
        % Calculate bivariate density only on first use and then cache it
        % so we can avoid this memory overhead.

        F2 = zeros(N, N);
        for q = 1:N
            for s = 1:N
                if q ~= s
                    d = dtmvnorm_marginal2([a(q), b(q), a(q), b(q)], [a(s), a(s), b(s), b(s)], q, s, zero_mean, sigma, lower, upper);
                    F2(q, s) = (d(1) - d(2)) - (d(3) - d(4));
                end
            end
        end
    

    % 2. Determine E[X_i, X_j]

    % Check if a[q] = -Inf or b[q]=+Inf, then F_a[q]=F_b[q]=0, but a[q] * F_a[q] = NaN and b[q] * F_b[q] = NaN
    F_a_q = ifelse(isinf(a), zeros(size(a)), a .* F_a); % n-dimensional vector q=1..N
    F_b_q = ifelse(isinf(b), zeros(size(b)), b .* F_b); % n-dimensional vector q=1..N

    for i = 1:N
        for j = 1:N
            sum = 0;
            for q = 1:N
                sum = sum + sigma(i, q) * sigma(j, q) * (sigma(q, q))^(-1) * (F_a_q(q) - F_b_q(q));
                if j ~= q
                    sum2 = 0;
                    for s = 1:N
                        % This term tt will be zero if the partial correlation coefficient \rho_{js.q} is zero!
                        % Even for s == q will the term be zero, so we do not need s!=q condition here
                        tt = (sigma(j, s) - sigma(q, s) * sigma(j, q) * (sigma(q, q))^(-1));
                        sum2 = sum2 + tt * F2(q, s);
                    end
                    sum2 = sigma(i, q) * sum2;
                    sum = sum + sum2;
                end
            end
            TVAR(i, j) = sigma(i, j) + sum;
            % General mean case: TVAR[i, j] = mean[j] * TMEAN[i] + mean[i] * TMEAN[j] - mean[i] * mean[j] + sigma[i, j] + sum
        end
    end

    % 3. Determine Variance Cov(X_i, X_j) = E[X_i, X_j] - E[X_i]*E[X_j] for (0, sigma)-case
    TVAR = TVAR - TMEAN * TMEAN';
    else
        TVAR = NaN;
    end

% 4. Shift back by +mean for (mu, sigma)-case
TMEAN = TMEAN + mean;

result = struct('tmean', TMEAN, 'tvar', TVAR);
end

%% helper function
function result = ifelse(condition, trueValue, falseValue)
    % Ensure that the inputs are valid
    if ~ismatrix(condition) || ~ismatrix(trueValue) || ~ismatrix(falseValue)
        error('All inputs must be matrices of the same size.');
    end
    
    % Apply the condition element-wise
    result = trueValue;
    result(~condition) = falseValue(~condition);
end