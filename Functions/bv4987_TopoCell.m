% bv4987_TopoCell.m
% author: Bryan Vu
% date: 11/10/19
clc; clear all; close all;

%% Physical parameters
Length = 500; % [m] aquifer length
Height = 50; % [m] aquifer thickness
K_hyd = 2e-6; % [m/s] hydraulic conductivity
Dh = 15; % [m] Regional slope of water table
x0 = Length/2; % [m] location of local maximum
dh = 3.2; % [m] height of local water table maximum
dw = 70; % [m] width of local maximum
Nh = 60; % [1] Number of head contours
Ns = 20; % [1] Number of streamlines
s = dw/4; % standard deviation of Gaussian maximum

%% Analytical Equation of h(x, Height) boundary
hb =@(x) Height + Dh*(1-x./Length) + dh*exp( -(x-x0).^2./(2*s^2) );

%% Build Grid
Grid.xmin = 0; Grid.xmax = Length; Grid.Nx = 500;
Grid.ymin = 0; Grid.ymax = Height; Grid.Ny = 150;
Grid = build_grid(Grid);

%% Define Boundary Conditions
Param.dof_dir = Grid.dof_ymax;
Param.dof_f_dir = Grid.dof_f_ymax;
Param.dof_neu = [];  % natural conditions
Param.dof_f_neu = [];

g = hb(Grid.xc)';  % heterogeneous dirichlet for ymax faces

%% Build Operators and Boundary
[D,G,I] = build_ops(Grid);
[B,N,fn] = build_bnd(Param, Grid, I);

%% Solve BVP
K = K_hyd*ones(Grid.Ny, Grid.Nx);
% isotropic and harmonic mean
Kd = comp_mean(K, 1, -1, Grid);
L = -D*Kd*G;

fs = spalloc(Grid.N, 1, 0);
h = solve_lbvp(L, fs+fn, B, g, N);

%% Compute Flux and Streamfunction
q = comp_flux(D, Kd, G, h, fs, Grid, Param);
[PSI, psi_min, psi_max] = comp_streamfun(q, Grid);

%% Stagnation Point
qx = reshape(q(1:Grid.Nfx), Grid.Ny, Grid.Nx+1);
qy = reshape(q(Grid.Nfx+1:end), Grid.Ny+1, Grid.Nx);

[Xc,Yc] = meshgrid(Grid.xc,Grid.yc);  % cell centers
[Xx,Yx] = meshgrid(Grid.xf,Grid.yc);  % x fluxes
[Xy,Yy] = meshgrid(Grid.xc,Grid.yf);  % y fluxes

qx_int = interp2(Xx,Yx,qx,Xc,Yc);
qy_int = interp2(Xy,Yy,qy,Xc,Yc);

f = abs(qx_int) + abs(qy_int);
% exclude boundaries
f([1,end],[1,end]) = max(f(:));  % fix!!!!
min_val = min(f(:));
[row, col] = find(f==min_val);

%% Plotting
figure(1); 
xx = linspace(0,Length,length(g))';
plot(xx, g, 'linewidth', 2)
title('Head on top boundary')

figure('Position', [250, 300, 1250, 350]); 

figure(2);
hold on;
contour(Xx, Yx, qx, [0 0], '-g', 'LineWidth', 2,'color',[.83,.26,.22])
contour(Xy, Yy, qy, [0 0], '-r', 'LineWidth', 2,'color',[.2,.6,1]) 
scatter3(Grid.xc(col), Grid.yc(row), f(row,col), 50, 'filled', 'k')
% Flownet 
h = reshape(h,[Grid.Ny,Grid.Nx]);
h_style = {'LineWidth'; 2};
psi_style = {'LineWidth'; 2};
plot_flownet(Nh, Ns, h, PSI, h_style, psi_style, Grid)

title('Nullclines')
legend('qx','qy','stag. pt.')
hold off;

% Flow Cells
figure('Position', [250, 300, 1250, 350]); 

[Xs,Ys] = meshgrid(Grid.xf,Grid.yf);  % cell corners for PSI
figure(3);
hold on;
plot_flownet(Nh, Ns, h, PSI, h_style, psi_style, Grid)
contour(Xs, Ys, PSI, [PSI(row,col), PSI(row,col)], '-r', 'LineWidth', 2)
scatter3(Grid.xc(col), Grid.yc(row), f(row,col), 50, 'filled', 'k')
hold off;
title('Flow Cells')

