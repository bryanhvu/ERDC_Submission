function [h] = solve_lbvp(L, f, B, g, N)
% author: Bryan Vu
% date: 9/23/2019
% Description
% Computes the solution $u$ to the linear differential problem given by
%
% $$\mathcal{L}(u)=f \quad x\in \Omega $$
%
% with boundary conditions
%
% $$\mathcal{B}(u)=g \quad x\in\partial\Omega$$.
%
% Input:
% L = matrix representing the discretized linear operator of size N by N,
% where N is the number of degrees of freedom
% f = column vector representing the discretized r.h.s. and contributions
% due to nonhomogeneous Neumann BC�s of size N by 1
% B = matrix representing the constraints arising from Dirichlet BC�s of
% size Nc by N
% g = column vector representing the non-homogeneous Dirichlet BC�s of size
% Nc by 1.
% N = matrix representing a orthonormal basis for the null-space of B and
% of size N by (N-Nc).
% Output:
% u = column vector of the solution of size N by 1
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
% >> Param.g = 0; % set bnd head
% >> [B,N,fn] = build_bnd(Param,Grid,I); % Build constraint matrix and
% % the basis for its nullspace
% >> L = -D*G; % Laplacian operator
% >> fs = spalloc(Grid.N,1,0); % r.h.s. (zero)
% >> h = solve_lbvp(L,fs+fn,B,Param.g,N); % Solve linear boundary value problem

% split h into particular and homogeneous solutions
hp = B' * ((B*B')\g);
ho = N * ((N'*L*N)\(N'*(f-L*hp)));
h = ho + hp;


