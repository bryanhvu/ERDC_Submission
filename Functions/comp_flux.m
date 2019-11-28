function [q] = comp_flux(D,Kd,G,h,fs,Grid,Param)
% author: Bryan Vu
% date: 09/30/2019
% Description:
% Computes the mass conservative fluxes across all boundaries from the
% residual of the compatability condition over the boundary cells.
% Note: Current implementation works for all cases where one face
% is assigned to each bnd cell. So corner cells must have
% natural BC’s on all but one face.
%
% Input:
% D = N by Nf discrete divergence matrix.
% Kd = Nf by Nf conductivity matrix.
% G = Nf by N discrete gradient matrix.
% h = N by 1 vector of flow potential in cell centers.
% fs = N by 1 right hand side vector containing only source terms.
% Grid = structure containing grid information.
% Param = structure containing problem parameters and information about BC’s
%
% Output:
% q = Nf by 1 vector of fluxes across all cell faces,
%
% Example call:
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
% >> Grid = build_grid(Grid);
% >> [D,G,I] = build_ops(Grid);
% >> L = -D*G; fs = ones(Grid.Nx,1);
% >> Param.dof_dir = Grid.dof_xmin;
% >> Param.dof_f_dir = Grid.dof_f_xmin;
% >> g = 0;
% >> Param.dof_neu = []; Param.dof_f_neu =[];
% >> [B,N,fn] = build_bnd(Param,Grid);
% >> h = solve_lbvp(L,fs+fn,B,g,N);
% >> q = comp_flux(D,1,G,h,fs,Grid,Param);

%% Compute interior fluxes
q = -Kd*(G*h);  % from Darcy's Law

%% Reconstruct boundary fluxes
% all non-natural conditions
dof_bnd = [Param.dof_dir Param.dof_neu]'; 
dof_f_bnd = [Param.dof_f_dir Param.dof_f_neu]'; 

% r = Lh - fs;
res_bnd = -D*Kd*G*h - fs; 

% check sign (direction) of flux on the boundary
sign = -ones(length(dof_bnd), 1);
subset = ismember(Grid.dof_f_xmin, dof_bnd); 
sign(Grid.dof_f_xmin(subset)) = 1;  % flux left to right chosen positive

q(dof_f_bnd) = sign.*res_bnd(dof_bnd).*Grid.V(dof_bnd)./Grid.A(dof_f_bnd);