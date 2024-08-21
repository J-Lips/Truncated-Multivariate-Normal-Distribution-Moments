function density = dtmvnorm_marginal2(xq, xr, q, r, mean, sigma, lower, upper, loga)
% Computation of the bivariate marginal density F_{q,r}(x_q, x_r) (q != r)
% of truncated multivariate normal distribution 
% following the works of Tallis (1961), Leppard and Tallis (1989)
%
% References:
% Tallis (1961): 
%   "The Moment Generating Function of the Truncated Multi-normal Distribution"
% Leppard and Tallis (1989): 
%   "Evaluation of the Mean and Covariance of the Truncated Multinormal"
% Manjunath B G and Stefan Wilhelm (2009): 
%   "Moments Calculation for the Doubly Truncated Multivariate Normal Distribution"
% input
    % xq mx1 query points on q dimension
    % xr mx1 query points on r dimension
    % q index for dimension q
    % r Index for Dimension r
    % mean nx1 
    % sigma nxn
    % lower nx1
    % upper nx1
    % loga take logarithm (default false)
    if nargin < 5
        mean = zeros(1, size(sigma, 1));
    end
    if nargin < 6
        sigma = eye(length(mean));
    end
    if nargin < 7
        lower = -Inf(1, length(mean));
    end
    if nargin < 8
        upper = Inf(1, length(mean));
    end
    if nargin < 9
        loga = false;
    end

    % dimensionality
    n = size(sigma, 1);

    % number of xq values delivered
    N = length(xq);

    % input checks
    if n < 2
        error('Dimension n must be >= 2!');
    end

    if length(mean) ~= size(sigma, 1)
        error('mean and sigma have non-conforming size');
    end

    if ~(ismember(q, 1:n) && ismember(r, 1:n))
        error('Indexes q and r must be integers in 1:n');
    end

    if q == r
        error('Index q must be different than r!');
    end

    % Skalierungsfaktor der gestutzten Dichte (Anteil nach Trunkierung)
    alpha = mvncdf(lower, upper, mean, sigma);

    if n == 2
        density = zeros(1, N);
        indOut = xq < lower(q) | xq > upper(q) | xr < lower(r) | xr > upper(r) | isinf(xq) | isinf(xr);
        density(~indOut) = dmvnorm([xq(~indOut); xr(~indOut)]', mean([q,r]), sigma([q,r],[q,r]), false) / alpha;
        if loga
            density = log(density);
        end
        return;
    end

    % standard deviation for normalisation
    SD = sqrt(diag(sigma));

    % normalised bounds
    lower_normalised = (lower - mean) ./ SD;
    upper_normalised = (upper - mean) ./ SD;

    xq_normalised = (xq - mean(q)) / SD(q);      % (N x 1)
    xr_normalised = (xr - mean(r)) / SD(r);      % (N x 1)

    % Computing correlation matrix R from sigma (matrix (n x n)): 
    % R = D % sigma * D with diagonal matrix D as sqrt(sigma)
    % same as cov2cor()
    D = diag(1 ./ sqrt(diag(sigma)));
    R = D * sigma * D;

    %
    % Determine (n-2) x (n-2) correlation matrix RQR
    %
    RQR = zeros(n-2, n-2);
    RINV = inv(R);
    WW = zeros(n-2, n-2);
    M1 = 0;
    for i = 1:n
        if i ~= q && i ~= r
            M1 = M1 + 1;
            M2 = 0;
            for j = 1:n
                if j ~= q && j ~= r
                    M2 = M2 + 1;
                    WW(M1, M2) = RINV(i,j);
                end
            end
        end
    end
    WW = inv(WW); % JOH WW = inv(WW(1:(n-2),1:(n-2)));
    for i = 1:(n-2)
        for j = 1:(n-2)
            RQR(i, j) = WW(i, j) / sqrt(WW(i,i) * WW(j,j));
        end
    end

    %
    % Determine bounds of integration vector AQR and BQR (n - 2) x 1
    %
    % lower and upper integration bounds
    AQR = zeros(N, n-2);                    
    BQR = zeros(N, n-2);
    M2 = 0;  % counter = 1..(n-2)
    for i = 1:n
        if i ~= q && i ~= r
            M2 = M2 + 1;
            BSQR = (R(q, i) - R(q, r) * R(r, i)) / (1 - R(q, r)^2);    
            BSRQ = (R(r, i) - R(q, r) * R(q, i)) / (1 - R(q, r)^2);    
            RSRQ = (1 - R(i, q)^2) * (1 - R(q, r)^2);
            RSRQ = (R(i, r) - R(i, q) * R(q, r)) / sqrt(RSRQ);         % partial correlation coefficient R[r,i] given q

            % lower integration bound
            AQR(:,M2) = (lower_normalised(i) - BSQR * xq_normalised - BSRQ * xr_normalised) / sqrt((1 - R(i, q)^2) * (1 - RSRQ^2));
            AQR(isnan(AQR(:,M2)), M2) = -Inf;

            % upper integration bound
            BQR(:,M2) = (upper_normalised(i) - BSQR * xq_normalised - BSRQ * xr_normalised) / sqrt((1 - R(i, q)^2) * (1 - RSRQ^2));
            BQR(isnan(BQR(:,M2)), M2) = Inf;
        end
    end

    % Correlation matrix for r and q
    R2 = [1, R(q,r); R(q,r), 1];

    sigma2 = sigma([q,r],[q,r]);            

    density = zeros(1, N);
    valid_indices = xq >= lower(q) & xq <= upper(q) & xr >= lower(r) & xr <= upper(r) & ~isinf(xq) & ~isinf(xr);
    
    prob = zeros(1, N);
    for i = 1:N
        if valid_indices(i)
            if (n - 2) == 1
                % univariate case
                prob(i) = normcdf(BQR(i,:), zeros(size(AQR(i,:))), RQR) - normcdf(AQR(i,:), zeros(size(AQR(i,:))), RQR);
            else
                prob(i) = mvncdf(AQR(i,:), BQR(i,:), zeros(size(AQR(i,:))), RQR);
            end
        end
    end
    
    density(valid_indices) = dmvnorm([xq(valid_indices); xr(valid_indices)]', mean([q,r]), sigma2);
    density(valid_indices) = density(valid_indices) .* prob(valid_indices) / alpha;

    if loga
        density = log(density);
    end
end