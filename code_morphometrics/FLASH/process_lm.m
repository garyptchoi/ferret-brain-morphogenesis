function [lm1a,lm2a,index] = process_lm(lm1,lm2)
% Process the landmarks.
%
% If you use this code in your own work, please cite the following paper:
% [1] P. T. Choi, K. C. Lam, and L. M. Lui, 
%     "FLASH: Fast Landmark Aligned Spherical Harmonic Parameterization for Genus-0 Closed Brain Surfaces."
%     SIAM Journal on Imaging Sciences, vol. 8, no. 1, pp. 67-94, 2015.
%
% Copyright (c) 2014-2022, Gary Pui-Tung Choi
% https://math.mit.edu/~ptchoi/

lm1a = [];lm2a = []; index = zeros(length(lm1a),1);

for i=1:length(lm1)
    [lm1i, lm2i] = landmark_modify(lm1{i}, lm2{i});
    lm1a = [lm1a,lm1i];
    lm2a = [lm2a,lm2i];
    index(i) = length(lm1a);
end

