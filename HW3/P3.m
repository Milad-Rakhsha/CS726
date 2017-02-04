clear
clc
x = struct('p', [-1.2, 5.2]);
trparams = struct('maxit',100,'delta',.01,'hatDelta',10,...
    'eta',.01,'Delta0',1,'toler',1.0e-6);


% Objective function
f = @(x) sin(pi*x(1)).*sin(pi*x(2));

% Gradient of Objective function
gradf = @(x)[pi.*cos(pi.*x(1)).*sin(pi.*x(2));pi.*cos(pi.*x(2)).*sin(pi.*x(1))];

% Hessian of Objective function
hessf = @(x)reshape([-pi.^2.*sin(pi.*x(1)).*sin(pi.*x(2)),pi.^2.*cos(pi.*x(1)).*cos(pi.*x(2)), ...
    pi.^2.*cos(pi.*x(1)).*cos(pi.*x(2)),-pi.^2.*sin(pi.*x(1)).*sin(pi.*x(2))],[2,2]);

DoglegTR(f,x,trparams)