% Compare real and simulated brains using spherical parameterization

addpath(genpath('code_morphometrics'));

%% load the real and simulated (stepwise) P16 brain data

load('data/P16_comparison_stepwise.mat');
% v1: the vertex coordinates of the real P16 brain surface
% f1: the triangulation of the real P16 brain surface
% lm1: vertex indices of labelled landmark curves on the real P16 brain
% bigtri1: the index of a triangle in f1 located at the central part
%         (for yielding a more symmetric spherical parameterization result)
% v2: the vertex coordinates of the simulated P16 brain surface
% f2: the triangulation of the simulated P16 brain surface
% lm2: vertex indices of labelled landmark curves on the simulated P16 brain
% bigtri2: the index of a triangle in f2 located at the central part
%         (for yielding a more symmetric spherical parameterization result)

% compute the curvature
curvature1 = tricurv(f1,v1);
mean_curv1 = curvature1.km;
Gauss_curv1 = curvature1.kg;

curvature2 = tricurv(f2,v2);
mean_curv2 = curvature2.km;
Gauss_curv2 = curvature2.kg;

% shape index
ShapeIndexHK1 = 2/pi*atan(mean_curv1./sqrt(mean_curv1.^2-Gauss_curv1));
ShapeIndexHK2 = 2/pi*atan(mean_curv2./sqrt(mean_curv2.^2-Gauss_curv2));


plot_mesh_lm(v1,f1,lm1,ShapeIndexHK1);
title('Real P16 brain')
view([-25 10])

plot_mesh_lm(v2,f2,lm2,ShapeIndexHK2);
title('Simulated P16 brain')
view([-25 10])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spherical parameterization 

% FLASH: fast landmark aligned spherical harmonic parameterization
% (Choi et al., SIAM J. Imaging Sci., 2015)
[map,s1,s2] = flash(v1,f1,lm1,v2,f2,lm2,10,bigtri1,bigtri2);

% The parameterization results
plot_mesh_lm(map,f1,lm1,threshold(mean_curv1,0.75)); 
title('Real P16')
view([140 0])

plot_mesh_lm(s2,f2,lm2,threshold(mean_curv2,0.75)); 
title('Simulated P16')
view([140 0])

%% Evaluate the similarity

F = TriScatteredInterp(s2(:,1),s2(:,2),s2(:,3),ShapeIndexHK2,'nearest');
interp_ShapeIndexHK2 = F(map(:,1),map(:,2),map(:,3));
I1 = ShapeIndexHK1;
I2 = interp_ShapeIndexHK2;

% similarity score
score = 1 - sum(abs(I1-I2))/(2*length(I1))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reconstruction using spherical harmonics (SH)

% maximum order for SH
max_order = 15; % 5 15 25

% SH reconstruction of the real brain
theta_brain_sph = acos(map(:,3)); 
phi_brain_sph = atan2(map(:,2),map(:,1)); 
phi_brain_sph(map(:,1)==0 & map(:,2)==0) = 0;
Y = Y_sph(length(map),max_order,theta_brain_sph,phi_brain_sph);
mult = (Y'*Y)\Y';
coeff = mult*v1;
v1_reconstruct = Y*coeff;

% SH reconstruction of the simulated brain
theta_brain_sph2 = acos(s2(:,3)); 
phi_brain_sph2 = atan2(s2(:,2),s2(:,1)); 
phi_brain_sph2(s2(:,1)==0 & s2(:,2)==0) = 0;
Y2 = Y_sph(length(s2),max_order,theta_brain_sph2,phi_brain_sph2);
mult2 = (Y2'*Y2)\Y2';
coeff2 = mult2*v2;
v2_reconstruct = Y2*coeff2;

%% plot the reconstruction results
plot_mesh(v1_reconstruct,f1);
view([-25 10]);

plot_mesh(v2_reconstruct,f2);
view([-25 10]);
