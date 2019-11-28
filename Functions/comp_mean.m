function [Kd] = comp_mean(K,kvkh,p,Grid)
% author: Bryan Vu
% date: 10/24/2019
% Description:
% Takes coefficient field, K, defined at the cell centers and computes the
% mean specified by the power, p and returns it in a sparse diagonal
% matrix, Kd.
%
% Input:
% K = Ny by Nx grid containing cell centered values
% p = power of the generalized mean
% 1 (arithmetic mean)
% -1 (harmonic mean)
% Grid = structure containing information about the grid.
% kvkh = ratio of vertical to horizontal conductivity (anisotropy)
% 
% Output:
% Kd = Nf by Nf diagonal matrix of power means at the cell faces.
%
% Example call:
% K = @(x) 1+x.^3;
% Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
% Grid = build_grid(Grid);
% Kd = comp_mean(K(Grid.xc),1,Grid);
% 

% just one value across, treat Kd as a constant
if (p == -1) || (p == 1)
  if (Grid.Nx == Grid.N) || (Grid.Ny == Grid.N) % 1D
     mean = zeros(Grid.N+1,1);
     mean(2:Grid.N) = 2*(K(1:Grid.N-1).^p + K(2:Grid.N).^p).^(1/p);
     Kd = spdiags(mean,0,Grid.N+1,Grid.N+1);
  elseif (Grid.N > Grid.Nx) || (Grid.N > Grid.Ny) % 2D
     % appropriating mean of Kh on cell centers to x faces
     mean_x = zeros(Grid.Ny,Grid.Nx+1);
     mean_x(:,2:Grid.Nx) = (.5*(K(:,1:Grid.Nx-1).^p + K(:,2:Grid.Nx).^p)).^(1/p);
     mean_x = reshape(mean_x,Grid.Nfx,1);
     
     % appropriating mean of Kv on cell centers to y faces
     K = K * kvkh; % account for anisotropy
     mean_y = zeros(Grid.Ny+1,Grid.Nx);
     mean_y(2:Grid.Ny,:) = (.5*(K(1:Grid.Ny-1,:).^p + K(2:Grid.Ny,:).^p)).^(1/p);
     mean_y = reshape(mean_y,Grid.Nfy,1);

     % Kd is Nf by Nf, so Nfx + Nfy on both sides
     mean = [mean_x; mean_y];
     Kd = spdiags(mean,0,Grid.Nf,Grid.Nf);
  else
    error('3D permeability is not implemented')
  end
else
  error('This power does not have significance.')
end

