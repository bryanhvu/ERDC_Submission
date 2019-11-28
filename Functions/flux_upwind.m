function [A] = flux_upwind(q,Grid)
% author: Bryan Vu
% date: 11/16/19
% Description:
% This function computes the upwind flux matrix from the flux vector.
%
% Input:
% q = Nf by 1 flux vector from the flow problem.
% Grid = structure containing all pertinent information about the grid.
%
% Output:
% A = Nf by N matrix containing the upwinded fluxes
% 
% Example call:
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
% >> Grid = build_grid(Grid);
% >> q = ones(Grid.Nf,1);

% >> [A] = flux_upwind(q,Grid);

if ((Grid.Nx>1) && (Grid.Ny==1)) || ((Grid.Nx==1) && (Grid.Ny>1)) % 1D
    %% 1D Advection Matrix
    qn = min(q(1:Grid.Nx),0);
    qp = max(q(2:Grid.Nx+1),0);
    
    A = spalloc(Grid.Nfx, Grid.Nx, 2*Grid.Nx);
    A(1:Grid.Nfx+1:end) = qn;
    A(2:Grid.Nfx+1:end) = qp;

elseif (Grid.Nx>1) && (Grid.Ny>1) % 2D
    % implement later
    error('2D not implemented')
else
    error('3D not implemented')
end