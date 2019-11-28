function [A] = flux_upwind(q,Grid)
% author: Bryan Vu
% date: November 15, 2019
% Description:
% This function computes the upwind flux matrix from the flux vector.
%
% Input:
% q = Nf by 1 flux vector from the flow problem.
% Grid = structure containing all pertinent information about the grid.
%
% Output:
% A = Nf by Nf matrix contining the upwinded fluxes
% 
% Example call:
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
% >> Grid = build_grid(Grid);
% >> q = ones(Grid.Nf,1);

% >> [A] = flux_upwind(q,Grid);
Nx = Grid.Nx; Ny = Grid.Ny; Nz = Grid.Nz; N = Grid.N;
Nfx = Grid.Nfx; % # of x faces
Nfy = Grid.Nfy; % # of y faces
if ((Nx>1) && (Ny==1)) || ((Nx==1) && (Ny>1)) % 1D
%% One dimensional
% code on notes
qn = min(q(1:Nx),0);
qp = max(q(2:Nx+1),0);
A = spalloc(Grid.Nfx, Grid.Nx, 2*Grid.Nx);
% A = zeros(Nfx, Nx);
A(1:Nfx+1:end) = qn(:);
A(2:Nfx+1:end) = qp(:);
% A = diag(qn) + diag(qp,-1);
elseif (Nx>1) && (Ny>1) % 2D
% Later
end