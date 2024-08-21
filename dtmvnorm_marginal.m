function density = dtmvnorm_marginal(xn, n, mean, sigma, lower, upper, loga)
% Density function for the marginal density of a truncated multivariate normal distribution
% Reference: Jack Cartinhour (1990) "One-dimensional marginal density functions of a truncated multivariate normal density function"
% input
    % xn Vector of length m of points at which the marginal density is evaluated
    % n Index (1..n) of the variable whose marginal density is to be computed
    % mean (nx1) Mean vector
    % sigma (nxn) Covariance matrix
    % lower, upper Truncation vectors lower <= x <= upper (nx1)
% output
    % mx1 marginal density at query points xn

    if nargin < 7
        loga = false;
    end
    if size(sigma, 1) ~= size(sigma, 2)
        error('sigma must be a square matrix');
    end

    if length(mean) ~= size(sigma, 1)
        error('mean and sigma have non-conforming size');
    end

    % Number of dimensions
    k = length(mean);

    if n < 1 || n > length(mean) || ~isnumeric(n) || length(n) > 1 || ~ismember(n, 1:length(mean))
        error('n must be an integer scalar in 1..length(mean)');
    end

    % Univariate case, see Greene (2003), p.573
    if k == 1
        prob = normcdf(upper, mean, sqrt(sigma)) - normcdf(lower, mean, sqrt(sigma));
        density = zeros(size(xn));
        density(lower <= xn & xn <= upper) = normpdf(xn(lower <= xn & xn <= upper), mean, sqrt(sigma)) / prob;
        if loga
            density = log(density);
        end
        return;
    end

    % Covariance matrix; after standardization, correlation matrix
    C = sigma;

    % Inverse covariance matrix, precision matrix
    A = inv(sigma);

    % Partitioning of A and C
    A_1 = A;
    A_1(:, n) = [];
    A_1(n, :) = [];
    A_1_inv = inv(A_1);

    C_1 = C;
    C_1(:, n) = [];
    C_1(n, :) = [];
    c_nn = C(n, n);
    c = C;
    c(n, :) = [];
    c = c(:, n);

    % Partitioning of mean vector mu
    mu = mean;
    mu_1 = mean;
    mu_1(n) = [];
    mu_n = mean(n);

    % Scaling factor of the density
    p = mvncdf(lower, upper, mu, C);

    f_xn = zeros(size(xn));
    for i = 1:length(xn)
        if (~(lower(n) <= xn(i) && xn(i) <= upper(n))) || isinf(xn(i))
            f_xn(i) = 0;
            continue;
        end

        % m(x_n) --> (n-1x1)
        % Be careful with m=c(Inf, Inf, NaN) and c=0
        m = mu_1 + (xn(i) - mu_n) * c / c_nn;

        % SW: Possibly optimize with vectorized version of pmvnorm() which accepts different bounds
        % for univariate density, pmvnorm() does not accept corr=
        f_xn(i) = exp(-0.5 * (xn(i) - mu_n)^2 / c_nn) * mvncdf(lower([1:n-1, n+1:end]), upper([1:n-1, n+1:end]), m, A_1_inv);
    end
    density = 1 / p / sqrt(2 * pi * c_nn) * f_xn;
    if loga
        density = log(density);
    end
end