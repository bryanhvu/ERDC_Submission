% Finding age of groundwater for 1D aquifer problem
clear variables; close all; clc;

%% add paths
addpath('../HW5')
addpath('../HW6')

%% Analytical Solution
% a_real =@(PHI, Q, x) (PHI./Q).*x; 
% syms a(x)
% a(x) = piecewise('x <= len', phi1.*x/q, 'x > len' phi2.*x/q)
% okay got a new function, but ode45 getting me crazy answer, look later

%% Parameters
yrs2sec = 3600*24*365;
km2m = 1e3;

Length = 20*km2m;
len = 8*km2m;
K1 = 5e-3;  % assign K from this
K2 = 5e-5; 
phi1 = 0.3;
phi2 = 0.1;

%% Grid
Grid.xmin = 0;
Grid.xmax = Length;
Grid.Nx = 100;

Grid = build_grid_2D(Grid);

%% Boundary Conditions
Param.dof_dir = [Grid.dof_xmax];
Param.dof_f_dir = [Grid.dof_f_xmax];
Param.dof_neu = [Grid.dof_xmin];
Param.dof_f_neu = [Grid.dof_f_xmin];

Param.qb = 1/yrs2sec;  % condition for neumann contribution, 1 m/s
g = 0;  % homogeneous

%% Build stuff
[D, G, I] = build_ops_2D(Grid);
[B, N, fn] = build_bnd(Param, Grid, I);

%% K
% K = zeros(Grid.Nx);
% K1p = (Grid.xc <= len) * K1;
% K2p = (Grid.xc > len) * K2;
% K = diag(K1p + K2p);   % don't mind rounding error, dumb

K = ones(Grid.Nx,1)*K1;
K(Grid.xc > len) = K2;

% isotropic and harmonic mean?
Kd = comp_mean_2D(K, 1, -1, Grid);
L1 = -D*Kd*G;

fs = spalloc(Grid.Nx, 1, 0);  % since no source term we're gucci
h = solve_lbvp(L1, fs+fn, B, g, N);
q = comp_flux(D, Kd, G, h, fs, Grid, Param);

%% advection
[A] = flux_upwind(q, Grid);
% linear operator
L2 = D*A;

% solve
phi = ones(Grid.Nx,1)*phi1;
phi(Grid.xc > len) = phi2;

%% ODE Solve
% dq = diff(q(:))./diff(Grid.xf(:));
% x = linspace(0, Length, 2*Grid.Nx);

% may have to input this into age_ode function
% q = interp1(Grid.xf, q, Grid.xc)';
% % dq = interp1(Grid.xc, dq, x)';
% % phi = interp1(Grid.xc, phi, x)';

% issue running ode solver
xspan = [0, Length];
ic = 0;
opts = odeset('RelTol', 1e-2, 'AbsTol', 1e-4);
sol = ode45(@(x, a) age_ode(x, a, q, phi, Grid), xspan, ic, opts);
% solution is NaN for about half of the problem, is q approaching 0?

% okay, soo taking a page from struct A and this class too
% potential solution involves a reduced system matrix
% where the singularities are removed. solve_lbvp already does that
% only thing now is to figure out proper BC on the age? or something

% update, with full grid looks good. Just need to fix end BC
% could be an easy fix. 

% a = L2 \ phi;
% g = ar(end) % analytical solution at x=L. 
a = solve_lbvp(L2, phi, B, g, N);

%% Plotting
% hold on;
plot(sol.x, sol.y/yrs2sec)
ylim([100
% xlim([0, 60])
% q = interp1(Grid.xf, q, Grid.xc, 'pchip');
% hold on;
% plot(Grid.xc/km2m, a_real(phi, q, Grid.xc))
% plot(Grid.xc/km2m, a/yrs2sec, 'linewidth', 2)
% hold off;
% xlabel('distance along aquifer: x[km]')
% ylabel('groundwater age: a [yrs]')