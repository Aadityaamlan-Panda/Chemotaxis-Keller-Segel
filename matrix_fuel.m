clear; clc;

% Parameters
Dc = 10e-5;
beta = 0.01;
U = 1/3 * 10^-4;
f0 = 0;
M0 = 1;
gamma = U * (1 - f0) / M0;

% Discretization
N = 100; % Number of grid points
r = linspace(0, 1, N)'; % Grid points
h = r(2) - r(1); % Step size

% Initialize unknowns
rho = zeros(N, 1);
S = zeros(N, 1);
f = zeros(N, 1);

% Boundary conditions
f(1) = f0; % f(0) = f0
f(end) = 1; % f(1) = 1

% Construct the finite difference matrix
A = zeros(3 * N, 3 * N); % System matrix
b = zeros(3 * N, 1); % Right-hand side

% Fill in the matrix A and vector b for the interior points
for i = 2:N-1
    % For rho equation
    A(i, i-1) = -1 / (2*h); % rho(i-1)
    A(i, i) = U + 1 / (2*h); % rho(i)
    A(i, N+i) = -rho(i); % S(i)
    A(i, i+1) = 1 / (2*h); % rho(i+1)
    
    % For S equation
    A(N+i, N+i-1) = -1 / (2*h); % S(i-1)
    A(N+i, N+i+1) = 1 / (2*h); % S(i+1)
    A(N+i, i) = -f(i); % rho(i)

    % For f equation
    A(2*N+i, 2*N+i-1) = -U / (2*h) + beta / (h^2); % f(i-1)
    A(2*N+i, 2*N+i) = -2 * beta / (h^2); % f(i)
    A(2*N+i, 2*N+i+1) = U / (2*h) + beta / (h^2); % f(i+1)
    A(2*N+i, i) = -gamma; % rho(i)
end

% Apply boundary conditions
A(1, 1) = 1; % rho(0) = 0
A(N, N) = 1; % rho(1) = 0

% Adjust for f boundary conditions
A(2*N+1, 2*N+1) = 1; % f(0) = f0
A(3*N, 3*N) = 1; % f(1) = 1

% Solve the linear system Ax = b
sol = A \ b;

% Extract solutions for rho, S, and f
rho = sol(1:N);
S = sol(N+1:2*N);
f = sol(2*N+1:end);

% Plot the results
figure;
plot(r, rho, '-r', 'DisplayName', '\rho(r)');
hold on;
plot(r, S, '-b', 'DisplayName', 'S(r)');
plot(r, f, '-g', 'DisplayName', 'f(r)');
xlabel('r');
ylabel('Solutions');
legend;
title('Finite Difference Solution for \rho(r), S(r), and f(r)');
