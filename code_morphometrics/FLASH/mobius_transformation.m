function vertex = mobius_transformation(vertex, p, q)

% Use three points q1, q2, q3 on sphere 2 to fix the transformation for sphere 1
% p,q: 3x3 matrix (each row = coordinates of 1 point)
% Trial:10.1.1, 9.1.1
% p1 = 40534 +1;q1 = 44481 +1;
% p2 = 34590 + 1;q2 = 10329 + 1;
% p3 = 8342 + 1; q3 = 8667 + 1;

numofvertex = size(vertex,1);

%% Project to complex plane
plane1 = complex(vertex(:,1)./(1-vertex(:,3)), vertex(:,2)./(1-vertex(:,3)));

point1 = complex(p(:,1)./(1-p(:,3)), p(:,2)./(1-p(:,3)));
point2 = complex(q(:,1)./(1-q(:,3)), q(:,2)./(1-q(:,3)));

p1 = point1(1,:);
p2 = point1(2,:);
p3 = point1(3,:);

q1 = point2(1,:);
q2 = point2(2,:);
q3 = point2(3,:);

%%
a = det([p1*q1,q1,1; p2*q2,q2,1; p3*q3,q3,1]);
b = det([p1*q1,p1,q1; p2*q2,p2,q2; p3*q3,p3,q3]);
c = det([p1,q1,1; p2,q2,1; p3,q3,1]);
d = det([p1*q1,p1,1; p2*q2,p2,1; p3*q3,p3,1]);

mobius_result = (a*plane1(:)+b)./(c*plane1(:)+d);
X = real(mobius_result(:));
Y = imag(mobius_result(:));

%% Map it to a sphere
vertex = [2*X./(1+X.^2+Y.^2) , 2*Y./(1+X.^2+Y.^2) , (-1+X.^2+Y.^2)./(1+X.^2+Y.^2)];

XX = [(1:numofvertex)',X];
XX_sort = sortrows(XX,2);
if isnan(XX_sort(numofvertex,2)) % the vertex is at (0,0,1) originally
    vertex(XX_sort(numofvertex,1),1:3) = [0.0001,0.0001,0.99999];
end

