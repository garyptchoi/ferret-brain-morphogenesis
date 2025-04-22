function vertex = mobius_transformation_ls(vertex, p, q)
% Least square method
% p, q are nx3 coorinates of landmarks

% If you use this code in your own work, please cite the following paper:
% [1] P. T. Choi, K. C. Lam, and L. M. Lui, 
%     "FLASH: Fast Landmark Aligned Spherical Harmonic Parameterization for Genus-0 Closed Brain Surfaces."
%     SIAM Journal on Imaging Sciences, vol. 8, no. 1, pp. 67-94, 2015.
%
% Copyright (c) 2013-2022, Gary Pui-Tung Choi
% https://math.mit.edu/~ptchoi/

%% Project to complex plane
plane = complex(vertex(:,1)./(1-vertex(:,3)), vertex(:,2)./(1-vertex(:,3)));

if size(p,2) == 3
    z = complex(p(:,1)./(1-p(:,3)), p(:,2)./(1-p(:,3)));
elseif size(p,2) == 2 
    z = complex(p(:,1), p(:,2));
else
    z = p;
end

if size(q,2) == 3
    w = complex(q(:,1)./(1-q(:,3)), q(:,2)./(1-q(:,3)));
elseif size(q,2) == 2
    w = complex(q(:,1), q(:,2));
else
    w = q;
end


g = 4./(1+abs(z).^2); 
x = [2*sum(g.*abs(z).^2), 0, sum(g.*(z+conj(z))), 1i*sum(g.*(conj(z)-z)); ...
    0, 2*sum(g.*abs(z).^2), 1i*sum(g.*(z-conj(z))), sum(g.*(conj(z)+z)); ...
    sum(g.*(z+conj(z))), 1i*sum(g.*(z-conj(z))), 2*sum(g), 0;...
    1i*sum(g.*(conj(z)-z)), sum(g.*(z+conj(z))), 0, 2*sum(g)]\...
    [sum(g.*(z.*conj(w)+conj(z).*w)); 1i*sum(g.*(z.*conj(w)-conj(z).*w)); sum(g.*(w+conj(w))); 1i*sum(g.*(conj(w)-w))];

a = x(1); b = x(2); c = x(3); d = x(4);
result = (a+1i*b)*plane + (c+1i*d);


vertex = [2*real(result)./(1+abs(result).^2), 2*imag(result)./(1+abs(result).^2), (abs(result).^2-1)./(1+abs(result).^2)];


