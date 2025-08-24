clear all
close all
clc

% Define constants
Dc = 10e-5;
beta = 0.01;
U = 1/3*10^-4;
f0 = 0;
M0 = 1;
gamma = U * (1 - f0) / M0;

% Define the initial guess for the solution
solinit = bvpinit(linspace(0, 1, 10), @guess);

% Solve the boundary value problem
sol = bvp4c(@odeFun, @bcFun, solinit);

% Extract solutions
r = linspace(0, 1, 100);
sol_rho = deval(sol, r, 1);
sol_drho = deval(sol, r, 2);
sol_c = deval(sol, r, 3);
sol_dc = deval(sol, r, 4);
sol_f = deval(sol, r, 5);
sol_df = deval(sol, r, 6);

% Plot solutions
figure;
plot(r, sol_rho, '-r', 'DisplayName', '\rho(r)');
hold on;
plot(r, sol_c, '-b', 'DisplayName', 'c(r)');
plot(r, sol_f, '-g', 'DisplayName', 'f(r)');
xlabel('r');
ylabel('Solutions');
legend;
title('Numerical Solutions for \rho(r), c(r), and f(r)');

% Function to define ODEs
function dydr = odeFun(r, y)
    global beta gamma U;
    dydr = zeros(6,1);
    
    % Assign variables for readability
    rho = y(1);
    drho = y(2);
    c = y(3);
    dc = y(4);
    f = y(5);
    df = y(6);
    
    % First-order system of ODEs
    dydr(1) = drho; % d(rho)/dr
    dydr(2) = -U * drho - rho * dc; % d2(rho)/dr2
    dydr(3) = dc; % d(c)/dr
    dydr(4) = -f * rho; % d2(c)/dr2
    dydr(5) = df; % d(f)/dr
    dydr(6) = (-U * df + beta * (df - df) - gamma * rho); % d2(f)/dr2 - fixed expression
end

% Function to define boundary conditions
function res = bcFun(ya, yb)
    global f0;
    res = [ya(1); yb(1); ya(3); yb(3)-1; ya(5)-f0; yb(5)-1];
end

% Initial guess for the solution
function yinit = guess(r)
    yinit = [0; 0; 0; 0; 0; 0];
end
