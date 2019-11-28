function [Grid] = build_grid(Grid)
% author: Bryan Vu bv4987
% date: 10/15/2019
% Description:
% This function takes in minimal definition of the computational domain 
% and computes all pertinent information about the grid (1D or 2D).
%
% Input: 
% Grid.xmin = left boundary of the domain
% Grid.xmax = right bondary of the domain
% Grid.Nx = number of grid cells
% Grid.ymin = lower boundary of the domain
% Grid.ymax = upper boundary of the domain
% Grid.Ny = number of grid cells in y direction
% Grid.psi_x0 = origin of streamfunction path integration
% Grid.psi_dir = direction of integration
%
% Output: 
% Grid.Lx = scalar length of the horizontal domain
% Grid.dx = scalar cell width
% Grid.Nfx = number of fluxes in x-direction
% Grid.xc = Nx by 1 column vector of cell center locations
% Grid.xf = Nfx by 1 column vector of cell face locations
% Grid.dof = Nx by 1 column vector from 1 to Nx containing the degrees of freedom, i.e. cell numbers
% Grid.dof_xmin = scalar cell degree of freedom corresponding to the left boundary
% Grid.dof_xmax = scalar cell degree of freedom corresponding to the right boundary
% Grid.dof_f_xmin = scalar face degree of freedom corresponding to the left boundary
% Grid.dof_f_xmax = scalar face degree of freedom corresponding to the right boundary
%
% Grid.Ly = scalar length of the vertical domain
% Grid.dy = scalar cell height
% Grid.Nfy = number of fluxes in y-direction
% Grid.yc = Ny by 1 column vector of cell center y locations
% Grid.yf = Nfy by 1 column vector of cell face y locations
% Grid.dof = Ny by 1 column vector from 1 to Ny containing the degrees of freedom, i.e. cell numbers
% Grid.dof_ymin = scalar cell degree of freedom corresponding to the left boundary
% Grid.dof_ymax = scalar cell degree of freedom corresponding to the right boundary
% Grid.dof_f_ymin = scalar face degree of freedom corresponding to the left boundary
% Grid.dof_f_ymax = scalar face degree of freedom corresponding to the right boundary
%
% Grid.dof_cell = scalar degree of freedom corresponding with cell value (2D)
% Grid.N = total cells
% Grid.Nf = total faces
% Grid.V = column vector of length N that contains cell volumes
% Grid.A = column vector of length Nf that contains all face areas
%
% Example call:
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
% >> Grid = build_grid(Grid);

if any(~isfield(Grid, {'Nx','Ny'}))  % 1D
    %% Grid Lengths and Cell Widths
    Grid.N = Grid.Nx;
    Grid.Lx = Grid.xmax - Grid.xmin;
    Grid.dx = Grid.Lx/Grid.Nx;
    Grid.Nfx = Grid.Nx + 1;  % # of faces for 1D case
    
    %% Cell and Face locations
    Grid.xc = (Grid.dx/2):Grid.dx:Grid.Lx;
    Grid.xf = 0:Grid.dx:Grid.Lx;

    %% Degrees of Freedom
    Grid.dof = 1:Grid.Nx;
    Grid.dof_f = 1:Grid.Nfx;
    Grid.dof_xmin = Grid.dof(1);
    Grid.dof_xmax = Grid.dof(end);
    Grid.dof_f_xmin = Grid.dof_f(1);
    Grid.dof_f_xmax = Grid.dof_f(end);

    %% Setting default values 
    Grid.Ny = 1; 
    Grid.Nfy = 1;
    Grid.Nz = 1;
    Grid.Nzy = 1;
    Grid.dy = 1; 
    Grid.dz = 1;
    
    %% Cell Volumes and Face Areas
    % equally spaced and sized cells
    Grid.V = Grid.dx.*Grid.dy.*Grid.dz.*ones(Grid.Nx,1);
    Grid.A = Grid.dy.*Grid.dz.*ones(Grid.Nfx,1);
    
elseif all(isfield(Grid, {'Nx', 'Ny'}))  % 2D
    %% Default Streamfunction Values
    % default location for origin is (xmin, ymin)
    if ~isfield(Grid,'psi_x0'); Grid.psi_x0 = 'xmin_ymin'; end
    % direction of integration starting x first then y
    if ~isfield(Grid,'psi_dir'); Grid.psi_dir = 'xy'; end

    %% Grid lengths and Cell widths
    Grid.Lx = Grid.xmax - Grid.xmin;
    Grid.dx = Grid.Lx/Grid.Nx;
    Grid.Ly = Grid.ymax - Grid.ymin;
    Grid.dy = Grid.Ly/Grid.Ny;

    %% Numbering 
    % rows of cells * column of faces 
    Grid.Nfx = Grid.Ny * (Grid.Nx + 1);
    % rows of faces * columns of cells
    Grid.Nfy = Grid.Nx * (Grid.Ny + 1);
    % calculate total faces and cells
    Grid.N = Grid.Nx * Grid.Ny;
    Grid.Nf = Grid.Nfx + Grid.Nfy;

    %% Cell and Face locations
    % cell centers and faces in different rows share same x-coord
    Grid.xc = (Grid.dx/2):Grid.dx:Grid.Lx;
    Grid.xf = 0:Grid.dx:Grid.Lx;
    Grid.yc = (Grid.dy/2):Grid.dy:Grid.Ly;
    Grid.yf = 0:Grid.dy:Grid.Ly;

    %% Degrees of Freedom
    Grid.dof_cell = reshape(1:Grid.N, Grid.Ny, Grid.Nx);
    Grid.dof_f_x = reshape(1:Grid.Nfx, Grid.Ny, Grid.Nx+1);
    Grid.dof_f_y = reshape(1:Grid.Nfy, Grid.Ny+1, Grid.Nx);
 
    Grid.dof_xmin = Grid.dof_cell(:,1);
    Grid.dof_xmax = Grid.dof_cell(:,end);
    Grid.dof_ymin = Grid.dof_cell(1,:);
    Grid.dof_ymax = Grid.dof_cell(end,:);

    Grid.dof_f_xmin = Grid.dof_f_x(:,1);
    Grid.dof_f_xmax = Grid.dof_f_x(:,end);
    % add x faces to simulate stacking
    Grid.dof_f_ymin = Grid.dof_f_y(1,:) + Grid.Nfx; 
    Grid.dof_f_ymax = Grid.dof_f_y(end,:) + Grid.Nfx;

    %% Cell Volumes and Face Areas
    Grid.dz = 1;  % default
    % equally spaced and sized cells
    Grid.V = Grid.dx.*Grid.dy.*Grid.dz.*ones(Grid.N,1);
    Grid.A = Grid.dy.*Grid.dz.*ones(Grid.Nf,1);
else
    error('Unknown Grid setup.')
end 
end