% test_ops2d_for.m
% author: Marc Hesse
% date: 26 Feb 2013, 11 Oct 2015
clc,  close all
clear Grid Param
% Problem parameters
nx = 1;
ny = 2;
h0 = 1;

% Define the computational grid
Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 50; 
Grid.ymin = 0; Grid.ymax = 1; Grid.Ny = 30; 
Grid = build_grid(Grid);


[Xc,Yc] = meshgrid(Grid.xc,Grid.yc);     % 2D coords of cell centers
[Xx,Yx] = meshgrid(Grid.xf,Grid.yc);     % 2D coords of x-faces
[Xy,Yy] = meshgrid(Grid.xc,Grid.yf);     % 2D coords of y-faces

% Generate operators
[D,G,I] = build_ops_2D(Grid);

%% Analytical derivatives;
h      = @(x,y)  cos(2*pi*nx/Grid.Lx*x).*cos(2*pi*ny/Grid.Ly*y);
dhdx   = @(x,y) -2*pi*nx/Grid.Lx*sin(2*pi*nx/Grid.Lx*x).*cos(2*pi*ny/Grid.Ly*y);
dhdy   = @(x,y) -2*pi*ny/Grid.Ly*cos(2*pi*nx/Grid.Lx*x).*sin(2*pi*ny/Grid.Ly*y);
d2hdx2 = @(x,y)  (2*pi*nx/Grid.Lx)^2*cos(2*pi*nx/Grid.Lx*x).*cos(2*pi*ny/Grid.Ly*y);
d2hdy2 = @(x,y)  (2*pi*ny/Grid.Ly)^2*cos(2*pi*nx/Grid.Lx*x).*cos(2*pi*ny/Grid.Ly*y);
Lap_h  = @(x,y)   d2hdx2(x,y) + d2hdy2(x,y);

%% Numerical derivatives
q = G*h(Xc(:),Yc(:));
QX = reshape(q(1:Grid.Nfx),Grid.Ny,Grid.Nx+1);
QY = reshape(q(Grid.Nfx+1:Grid.Nf),Grid.Ny+1,Grid.Nx);
Lap_num = reshape(-D*q,Grid.Ny,Grid.Nx);

figure('name','Testing 2D Operators')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot 331
contourf(Xx,Yx,dhdx(Xx,Yx)), colorbar
axis([0 Grid.Lx 0 Grid.Ly])
xlabel 'x', ylabel 'y'
axis equal tight
title 'dh/dx'

subplot 332
contourf(Xx,Yx,QX), colorbar
axis([0 Grid.Lx 0 Grid.Ly])
xlabel 'x', ylabel 'y'
axis equal tight
title 'qx'

subplot 333
temp = dhdx(Xx,Yx);
err = (dhdx(Xx,Yx) - QX)/max(temp(:));
contourf(Xx,Yx,err), colorbar
axis([0 Grid.Lx 0 Grid.Ly])
xlabel 'x', ylabel 'y'
axis equal tight
title 'Error x'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot 334
contourf(Xy,Yy,dhdy(Xy,Yy)), colorbar
axis([0 Grid.Lx 0 Grid.Ly])
xlabel 'x', ylabel 'y'
axis equal tight
title 'dh/dy'

subplot 335
contourf(Xy,Yy,QY), colorbar
axis([0 Grid.Lx 0 Grid.Ly])
xlabel 'x', ylabel 'y'
axis equal tight
title 'qy'

subplot 336
temp = dhdy(Xy,Yy);
err = (temp - QY)/max(temp(:));
contourf(Xy,Yy,err), colorbar
axis([0 Grid.Lx 0 Grid.Ly])
xlabel 'x', ylabel 'y'
axis equal tight
title 'Error y'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot 337
contourf(Xc,Yc,Lap_h(Xc,Yc)), colorbar
axis([0 Grid.Lx 0 Grid.Ly])
xlabel 'x', ylabel 'y'
axis equal tight
title '\nabla^2h'

subplot 338
contourf(Xc,Yc,Lap_num), colorbar
axis([0 Grid.Lx 0 Grid.Ly])
xlabel 'x', ylabel 'y'
axis equal tight
title 'L*h'

subplot 339
temp = Lap_h(Xc,Yc);
err = (temp - Lap_num)/max(temp(:));
contourf(Xc,Yc,err), colorbar
axis([0 Grid.Lx 0 Grid.Ly])
xlabel 'x', ylabel 'y'
axis equal tight
title 'Error y'
