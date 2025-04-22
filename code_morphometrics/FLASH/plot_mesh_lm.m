function plot_mesh_lm(v,f,lm,quantity)

% Plot a mesh with landmarks.
% 
% Input: 
% v: nv x 3 vertex coordinates
% f: nf x 3 triangulations
% lm: nl x 1 cell array of landmark indices
% (optional) quantity: nv x 1 quantity defined on vertices
%
% If you use this code in your own work, please cite the following paper:
% [1] P. T. Choi, K. C. Lam, and L. M. Lui, 
%     "FLASH: Fast Landmark Aligned Spherical Harmonic Parameterization for Genus-0 Closed Brain Surfaces."
%     SIAM Journal on Imaging Sciences, vol. 8, no. 1, pp. 67-94, 2015.
%
% Copyright (c) 2014-2022, Gary Pui-Tung Choi
% https://math.mit.edu/~ptchoi/

if nargin == 4
    plot_mesh(v,f,quantity);
    colormap('Copper');
    shading interp;
    colorbar('off');
else
    plot_mesh(v,f);
end
axis('off');
set(gcf,'color','w'); 
view([90,10]);
hold on; 
if size(v,2) == 3
    for i = 1:length(lm)
      plot3(v(lm{i},1),v(lm{i},2),v(lm{i},3),'r-', 'LineWidth',5);
      plot3(v(lm{i},1),v(lm{i},2),v(lm{i},3),'o','MarkerSize',5,...
          'MarkerEdgeColor','r','MarkerFaceColor',[1,0,0]);
    end
else
    for i = 1:length(lm)
      plot(v(lm{i},1),v(lm{i},2),'r-', 'LineWidth',5);
      plot(v(lm{i},1),v(lm{i},2),'o','MarkerSize',5,...
          'MarkerEdgeColor','r','MarkerFaceColor',[1,0,0]);
    end
end
    