function [D,G,I]=build_ops(Grid)
% author: Bryan Vu bv4987
% date: 10/16/2019
% description:
% This function computes the discrete divergence and gradient matrices on a
% regular staggered grid using central difference approximations. The
% discrete gradient assumes homogeneous boundary conditions.
% Input:
% Grid = structure containing all pertinent information about the grid.
% Output:
% D = N by Nf discrete divergence matrix
% G = Nf by N discrete gradient matrix
% I = N by N identity matrix
%
% Example call:
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
% >> Grid = build_grid(Grid);
% >> [D,G,I]=build_ops(Grid);

I = speye(Grid.N);  % built the same for 1D and 2D
if any(~isfield(Grid,{'Nx', 'Ny'}))  % 1D
    %% Build 1D Operators 
    diag = (1/Grid.dx) * ones(Grid.N, 1);
    D = spdiags([-diag, diag], 0:1, Grid.N, Grid.Nfx);   
    G = -D';  % adjoints
    G(1,1) = 0; G(end,end) = 0;  % Using natural BC's
    
elseif all(isfield(Grid, {'Nx', 'Ny'}))  % 2D
    %% Assemble Components of Divergence Operator
    % Divergence Operator in y cells
    diag = (1/Grid.dy) * ones(Grid.Ny, 1);
    Dy1 = spdiags([-diag, diag], 0:1, Grid.Ny, Grid.Ny+1);
    Ix = speye(Grid.Nx);
    Dy = kron(Ix, Dy1);

    % Divergence Operator in x cells
    diag = (1/Grid.dx) * ones(Grid.Nx, 1);
    Dx1 = spdiags([-diag, diag], 0:1, Grid.Nx, Grid.Nx+1);
    Iy = speye(Grid.Ny);
    Dx = kron(Dx1, Iy);
    
    %% Construct 2D Operators
    D = [Dx, Dy];
    G = -D';  % adjoints

    % Using natural BC's
    G(Grid.dof_f_xmin,Grid.dof_xmin) = 0;
    G(Grid.dof_f_xmax,Grid.dof_xmax) = 0;
    G(Grid.dof_f_ymin,Grid.dof_ymin) = 0;
    G(Grid.dof_f_ymax,Grid.dof_ymax) = 0;
else
    error('Unknown Grid setup')
end
end
