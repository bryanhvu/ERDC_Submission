function [B,N,fn] = build_bnd(Param,Grid,I)
% author: Bryan Vu
% date: 9/20/2019
% Description:
% This function computes the operators and r.h.s vectors for both Dirichlet
% and Neumann boundary conditions.
%
% Input:
% Grid = structure containing all pertinent information about the grid.
% Param = structure containing all information about the physical problem
% in particular this function needs the fields:
% Param.dof_dir = Nc by 1 column vector containing
% the dof’s of the Dirichlet boundary.
% Param.dof_neu = column vector containing
% the dof’s of the Neumann boundary.
% Param.qb = column vector of prescribed fluxes on Neuman bnd.
% I = identity matrix in the full solution space
%
% Output:
% B = Nc by N matrix of the Dirichlet constraints
% N = (N-Nc) by (N-Nc) matrix of the nullspace of B
% fn = N by 1 r.h.s. vector of Neumann contributions
%
% Example call:
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
% >> Grid = build_grid(Grid);
% >> [D,G,I]=build_ops(Grid);
% >> Param.dof_dir = Grid.dof_xmin; % identify cells on Dirichlet bnd
% >> Param.dof_f_dir = Grid.dof_f_xmin; % identify faces on Dirichlet bnd
% >> Param.dof_neu = Grid.dof_xmax; % identify cells on Neumann bnd
% >> Param.dof_f_neu = Grid.dof_f_xmax; % identify cells on Neumann bnd
% >> Param.qb = 1; % set bnd flux
% >> [B,N,fn] = build_bnd(Param,Grid,I);

if ~isfield(Param, 'qb')
    Param.qb = 0;  % default
end

%% Compute Fluxes
fn = spalloc(Grid.N,1,0);
fn(Param.dof_neu) = Param.qb*Grid.A(Param.dof_f_neu)/Grid.V(Param.dof_neu);
% fn(Param.dof_neu) = Param.qb*Grid.V(Param.dof_f_neu)/Grid.A(Param.dof_neu);
% figure out which is which already 

%% Construct Constraint Space and Projection Matrix
dof_bnd = Param.dof_dir;
B = I(dof_bnd,:);
N = I;
N(:,dof_bnd) = [];