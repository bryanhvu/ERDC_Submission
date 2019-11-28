% bv4987_GWage1D.m
% author: Bryan Vu
% date: 11/18/19
clear variables; close all; clc;

%% Add Path to External Functions
addpath('../Functions')

%% Parameters
yrs2sec = 3600*24*365;
km2m = 1e3;

Length = 20*km2m;  % [m] length of aquifer
len = 8*km2m;  % [m] distance where piecewise changes occur
K1 = 5e-3;  % [1] hydraulic conductivity of the first segment
K2 = 5e-5;  % [1] hydraulic conductivity of the second segment
phi1 = 0.3;  % [1] porosity of the first segment
phi2 = 0.1;  % [1] porosity of the second segment

%% Build Grid
Grid.xmin = 0; 
Grid.xmax = Length; 
Grid.Nx = 100;
Grid = build_grid(Grid);

%% Define Boundary Conditions
Param.dof_dir = [Grid.dof_xmax];
Param.dof_f_dir = [Grid.dof_f_xmax];
Param.dof_neu = [Grid.dof_xmin];
Param.dof_f_neu = [Grid.dof_f_xmin];

Param.qb = 1/yrs2sec; 
g = 0;  % homogeneous

%% Build Operators and Boundary
[D, G, I] = build_ops(Grid);
[B, N, fn] = build_bnd(Param, Grid, I);

% variable hydraulic conductivity
K = ones(Grid.Nx,1)*K1;
K(Grid.xc > len) = K2;
% isotropic and harmonic mean
Kd = comp_mean(K, 1, -1, Grid);
L1 = -D*Kd*G;

%% Solve Flow Problem
fs = spalloc(Grid.N, 1, 0);
h = solve_lbvp(L1, fs+fn, B, g, N);
q = comp_flux(D, Kd, G, h, fs, Grid, Param);

%% Solve Analytical Advection Problem
% variable porosity
phi = ones(Grid.Nx, 1)*phi1;
phi(Grid.xc > len) = phi2;

xspan = [0, Length];
ic = 0;
opts = odeset('RelTol', 1e-2, 'AbsTol', 1e-4);
sol = ode45(@(x, a) age_ode(x, a, q, phi, Grid), xspan, ic, opts);

%% Solve Numerical Advection Problem
A = flux_upwind(q, Grid);
% linear operator for advection problem
L2 = D*A;

g = interp1(sol.x, sol.y, Grid.xc(end)); % analytical solution at x=L. 
a = solve_lbvp(L2, phi, B, g, N);

%% Plotting
hold on;
plot(sol.x/km2m, sol.y/yrs2sec, 'linewidth', 1.5)
plot(Grid.xc/km2m, a'/yrs2sec, '--r', 'linewidth', 1.5)
hold off;

xlabel('distance along aquifer: x[km]')
ylabel('groundwater age: a [yrs]')
title('Groundwater Age vs. Distance along Aquifer')
legend('Analytical', 'Numerical')