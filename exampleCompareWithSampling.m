% Example: comparison between numerical and analytical moment calculation
% numerical calculation uses the Truncated Multivariate Normal
% Generator by Zdravko Botev
% (https://www.mathworks.com/matlabcentral/fileexchange/53792-truncated-multivariate-normal-generator)

% in this example, a multidimensional normal distribution is constrained in
% both directions for all variables, the new mean and std are calculated
% and plotted. This is done numerically (through sampling of the
% truncated distribution) and analytically (using formulas from this
% package)
%%
% == problem definition ==
rng('default')

n = 8;
mu = rand(n,1); %mean
[i, j] = meshgrid(1:n, 1:n);
sigma = 1 * exp(-1 * (i - j).^2); %exponentially decaying corr-matrix

lower = mu - rand(n,1);
upper = mu + rand(n,1);

% == analytical calculation ==
resultFormula = mtmvnorm(mu, sigma, lower, upper, true);
% assign 1D std for each variable
resultFormula.tstd = real(sqrt(diag(resultFormula.tvar))); 

% == numerical calculation ==
% sample mean-free truncated multivariate normal and add mean
x = mvrandn(lower-mu, upper-mu,sigma,numRealizations)+mu;
% calculate mean and std from data
resultSampling.tmean = mean(x, 2);
resultSampling.tstd = std(x, [], 2);

% == visualization ==
figure();
hold on;
plot(mu, 'Color', [0.2 0.2 0.2], 'DisplayName', 'untrunc. mean')
plot(lower, 'Color', [0.7 0.7 0.7], 'DisplayName', 'lower/upper')
plot(upper, 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off')

result = resultSampling;
plot(result.tmean, 'Color', [1 0 0], 'LineStyle', '-', 'DisplayName', 'trunc. mean (sampling)')
plot(result.tmean + result.tstd, 'Color', [1 0.7 0.7], 'LineStyle', '-', 'DisplayName', 'trunc. mean +- std (sampling)')
plot(result.tmean - result.tstd, 'Color', [1 0.7 0.7], 'LineStyle', '-', 'HandleVisibility', 'off')

result = resultFormula;
plot(result.tmean, 'Color', [0 0 1], 'LineStyle', '--', 'DisplayName', 'trunc. mean (analytical)')
plot(result.tmean + result.tstd, 'Color', [0.7 0.7 1], 'LineStyle', '--', 'DisplayName', 'trunc. mean +- std (analytical)')
plot(result.tmean - result.tstd, 'Color', [0.7 0.7 1], 'LineStyle', '--', 'HandleVisibility', 'off')
legend('Location','northoutside')