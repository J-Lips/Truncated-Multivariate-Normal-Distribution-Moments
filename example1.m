% This is the example presented as Numerical Example 1 in  
% Moments Calculation For the Doubly Truncated Multivariate Normal Density
% B.G. Manjunath, S. Wilhelm 2021
% DOI: 10.35566/jbds/v1n1/p2
%% 
rng('default')

% == problem definition ==
mu = [0.5 0.5]';
sigma = [1  1.2; 
        1.2 2];
a = [-1 -inf]';
b = [0.5 1]';

% == calculation of moments ==
moments = mtmvnorm(mu, sigma, a, b, true);
disp('mean = ');disp(moments.tmean);
disp('Sigma = ');disp(moments.tvar);

% == marginal distribution ==
xquery1 = a(1)-0.1:0.01:b(1)+0.1;
F_1 = dtmvnorm_marginal(xquery1, 1, mu, sigma, a, b);
xquery2 = -3.5:0.01:b(2)+0.1;
F_2 = dtmvnorm_marginal(xquery2, 2, mu, sigma, a, b);

% == bivariate marginal distribution ==
[x1Grid, x2Grid] = meshgrid(xquery1, xquery2);
F12 = dtmvnorm_marginal2(reshape(x1Grid, 1, []), reshape(x2Grid, 1, []), 1, 2, mu, sigma, a, b);
F12 = reshape(F12, size(x1Grid));
%% == visualization ==
% === scatter plot with truncation ===
% generate points from the untruncated distribution
numRealizations = 3000;
x = mvnrnd(mu,sigma,numRealizations)';

% Plot the scatter plot of the points
figure('Position',[100 500 600 250]);
title('bivariate marginal density (x_1,x_2)');
hold on;
mask = true(size(x(1,:)));
for i = 1:size(x,1)
    mask = mask & (x(i,:)>=a(i)) &  (x(i,:)<=b(i));
end
scatter(x(1, mask), x(2, mask), 8, 'MarkerEdgeColor',[0 0 0]); 
scatter(x(1, ~mask), x(2, ~mask), 8, 'MarkerEdgeColor',0.7*[1 1 1]); 
x_rect = [a(1), b(1), b(1), a(1), a(1)];
y_rect = [a(2), a(2), b(2), b(2), a(2)];
x_rect(x_rect==-inf) = -5;
y_rect(y_rect==-inf) = -5;
plot(x_rect, y_rect, 'r--', 'LineWidth', 1);
xlabel('x1');
xlim([-4.5 4.5]);
ylabel('x2');
ylim([-4.5 4.5]);

contour(x1Grid, x2Grid, F12, 'LineWidth', 1.5, 'ShowText','on', 'EdgeColor','r', 'EdgeAlpha',0.8,  'TextStep',0.02);
hold off;

% === marginal density plots ===
figure('Position', [100 100 250 300]);
title('Marginal density x_1');
hold on;
xi = x(1, mask);
[f, xi] = ksdensity(xi);
plot(xi(xi>=a(1)&xi<=b(1)), f(xi>=a(1)&xi<=b(1)), 'k-');
plot(xquery1, F_1, 'r-');
%plot([a(1) a(1)], [0 1], 'r--');
%plot([b(1) b(1)], [0 1], 'r--');
ylabel('Density')
xlabel('x_1')
xlim([xquery1(1) xquery1(end)]);
ylim([0.15 0.95])

figure('Position', [400 100 250 300]);
title('Marginal density x_2');
hold on;
xi = x(2, mask);
[f, xi] = ksdensity(xi);
plot(xi(xi>=a(2)&xi<=b(2)), f(xi>=a(2)&xi<=b(2)), 'k-');
plot(xquery2, F_2, 'r-');
%plot([a(1) a(1)], [0 1], 'r--');
%plot([b(1) b(1)], [0 1], 'r--');
ylabel('Density')
xlabel('x_2')
xlim([xquery2(1) xquery2(end)]);
ylim([-0.05 0.55])
