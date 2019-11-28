function [] = plot_flownet(Nh, Ns, h, PSI, h_style, psi_style, Grid)
% author: Bryan Vu
% date: 10/18/2019
% Description: Plots a flownet with Nh equally spaced head contours
% and Ns equally spaced streamlines.
% 
% Input: 
% Nh = number of head contours
% Ns = number of streamlines
% h = N by 1 column vector of heads
% PSI = Ny+1 by Nx+1 matrix containing the stream function
% h_style = cell array specifying the linestyle for the head contours
% psi_style = cell array specifying the linestyle for the streamlines
% Grid = structure containing all information about the grid.
% 
% Output: none

[Xc,Zc] = meshgrid(Grid.xc, Grid.yc);
[Xf,Zf] = meshgrid(Grid.xf, Grid.yf);
h = reshape(h,Grid.Ny, Grid.Nx);

contour(Xf, Zf, PSI, Ns, psi_style{:})  % streamlines
contour(Xc, Zc, h, Nh, h_style{:})  % head contours

xlabel('x (m)'); 
ylabel('z (m)');