clear all
close all
clc
syms rho(r) f(r) S(r) beta Dc U f0 M0 gamma
Dc = 10e-5;
beta = 0.01;
U = 1/3*10^-4;
f0 = 0;
M0 = 1;
gamma = U*(1-f0)/M0;
eqns = [-U*rho == diff(rho,r)-rho*S, diff(S,r) == -f*rho, -U*diff(f,r) == -gamma*rho+beta*diff(f,r,2)];
conds = [rho(0) == 0, rho(1) == 0, f(0) == f0, f(1) == 1];
sol = dsolve(eqns,conds);