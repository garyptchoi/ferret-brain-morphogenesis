function [newedge1, newedge2] = landmark_modify(edge1,edge2)
% Modify two landmarks with different lengths.
% edge1, edge2: two row vectors 
%
% If you use this code in your own work, please cite the following paper:
% [1] P. T. Choi, K. C. Lam, and L. M. Lui, 
%     "FLASH: Fast Landmark Aligned Spherical Harmonic Parameterization for Genus-0 Closed Brain Surfaces."
%     SIAM Journal on Imaging Sciences, vol. 8, no. 1, pp. 67-94, 2015.
%
% Copyright (c) 2013-2022, Gary Pui-Tung Choi
% https://math.mit.edu/~ptchoi/

L1 = length(edge1);
L2 = length(edge2);

if L1 < L2
    seq = zeros(1,L1);
    for i = 1 : L1
        seq(i) = floor(i*L2/L1);
    end
    newedge1 = edge1;
    newedge2 = edge2(seq);

elseif L1 > L2
    seq = zeros(1,L2);
    for i = 1 : L2
        seq(i) = floor(i*L1/L2);
    end
    newedge1 = edge1(seq);
    newedge2 = edge2;
    
else
    newedge1 = edge1;
    newedge2 = edge2;
end