% This is the example presented as Numerical Example 2 in  
% Moments Calculation For the Doubly Truncated Multivariate Normal Density
% B.G. Manjunath, S. Wilhelm 2021
% DOI: 10.35566/jbds/v1n1/p2
%%
rng('default')

% == problem definition ==
mu = [0 0 0]';
sigma = [1.1   1.2   0; 
         1.2    2  -0.8;
          0   -0.8   3];
a = [-1 -inf -inf]';
b = [0.5 inf inf]';

% == calculation of moments ==
moments = mtmvnorm(mu, sigma, a, b, true);
disp('mean = ');disp(moments.tmean);
disp('Sigma = ');disp(moments.tvar);