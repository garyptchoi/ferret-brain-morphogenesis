function [map,s1,s2,f1_filled,f2_filled] = flash(v1,f1,lm1,v2,f2,lm2,lambda,bigtri1,bigtri2)
% Main program for the landmark-aligned spherical harmonic map for brain
% surfaces, including the initial sulcal landmark curve processing and the
% mapping computation. 
% Use flash.m instead if there is no need to process the landmark curves.
%
% Input:
% v1: mx3 matrix for vertices of the source mesh
% f1: px3 matrix for faces of the source mesh
% lm1: kx1 cell array with each entry containing a row vector of indices of the landmark curves on the source mesh
% v2: nx3 matrix for vertices of the target mesh
% f2: qx3 matrix for vertices of the target mesh
% lm2: kx1 cell array with each entry containing a row vector of indices of the landmark curves on the target mesh
% lambda: a non-negative landmark matching factor 
% (Optional) bigtri1: ID of the face to be "punctured" for spherical map of the source surface
% (Optional) bigtri2: ID of the face to be "punctured" for spherical map of the target surface
% If bigtri1 and bigtri 2 are not specified, the algorithm will look for the most regular triangle on each mesh. 
% Choosing a triangle too close to landmark curves may affect the computation. 
%
% Output:
% map: mx3 matrix for vertices of the landmark aligned spherical harmonic map
% s1: mx3 matrix for vertices of the spherical conformal parameterization of the source mesh
% s2: nx3 matrix for vertices of the spherical conformal parameterization of the target mesh
% f1_filled: px3 matrix for faces of the first mesh, 
%            possibly with holes filled to make it genus-0
% f2_filled: qx3 matrix for faces of the second mesh, 
%            possibly with holes filled to make it genus-0
%
% If you use this code in your own work, please cite the following paper:
% [1] P. T. Choi, K. C. Lam, and L. M. Lui, 
%     "FLASH: Fast Landmark Aligned Spherical Harmonic Parameterization for Genus-0 Closed Brain Surfaces."
%     SIAM Journal on Imaging Sciences, vol. 8, no. 1, pp. 67-94, 2015.
%
% Copyright (c) 2013-2022, Gary Pui-Tung Choi
% https://math.mit.edu/~ptchoi/

if length(lm1)~=length(lm2) || isempty(lm1) || isempty(lm2)
    error('The number of landmark curves are different.');
end

if lambda < 0
    error('lambda should be a nonnegative number.');
end

bd1 = meshboundaries(f1);
bd2 = meshboundaries(f2);

if ~isempty(bd1)
    warning(['Mesh 1 contains ', ...
        num2str(length(bd1)), ...
        ' boundaries instead of 0. The extra holes will be filled.']);
    f1 = hole_filling(f1,0);
end

if ~isempty(bd2)
    warning(['Mesh 2 contains ', ...
        num2str(length(bd2)), ...
        ' boundaries instead of 0. The extra holes will be filled.']);
    f2 = hole_filling(f2,0);
end

f1_filled = f1;
f2_filled = f2;

% to avoid error in mobius transformation
if length(lm1) == 1
    lm1_temp = cell2mat(lm1);
    lm2_temp = cell2mat(lm2);
    
    lm1 = cell(2,1); 
    lm1{1} = lm1_temp(1:length(lm1_temp)/2);
    lm1{2} = lm1_temp(round((length(lm1_temp)+1)/2):end);
    
    lm2 = cell(2,1); 
    lm2{1} = lm2_temp(1:length(lm2_temp)/2);
    lm2{2} = lm2_temp(round((length(lm2_temp)+1)/2):end);
end

% ensure equal number of landmark points on each curve
[lm1a,lm2a,index] = process_lm(lm1,lm2);

v0 = v1;
fprintf('Landmark resampling completed.\n');

%% obtain the conformal parameterizations
fprintf('Computing the spherical conformal map of the source mesh... ');
if exist('bigtri1','var')
    [s1,bigtri1] = spherical_conformal_map(v1,f1,bigtri1); 
else
    [s1,bigtri1] = spherical_conformal_map(v1,f1);
end
fprintf('Done.\n');

fprintf('Computing the spherical conformal map of the target mesh... ');
if exist('bigtri2','var')
    [s2,~] = spherical_conformal_map(v2,f2,bigtri2); 
else
    [s2,~] = spherical_conformal_map(v2,f2);
end
fprintf('Done.\n');

fprintf('Initial landmark mismatch energy = %d\n', ...
    sum((s2(lm2a,1)-s1(lm1a,1)).^2+(s2(lm2a,2)-s1(lm1a,2)).^2 + (s2(lm2a,3)-s1(lm1a,3)).^2));

%% Landmark-constrained optimized harmonic map on the plane
fprintf('Computing the landmark-constrained optimized harmonic map... ');
DistfromNP = 1-s1(lm1a(unique(sort([1,index,index(length(index)-1)+1]))),3);
DistfromNP = sortrows([(1:length(DistfromNP))',DistfromNP],2); 
lm1_endpt = lm1a(unique(sort([1,index,index(length(index)-1)+1])));
lm2_endpt = lm2a(unique(sort([1,index,index(length(index)-1)+1])));
% Apply Mobius transformation to roughly align the landmarks
sm1 = mobius_transformation_ls(s1, s1(lm1a,:), s2(lm2a,:));

numofvertex1 = size(sm1, 1);
[~,NorthTri1] = min(abs(((sm1(f1(:,1),1)+sm1(f1(:,2),1)+sm1(f1(:,3),1)))/3) + abs(((sm1(f1(:,1),2)+sm1(f1(:,2),2)+sm1(f1(:,3),2)))/3) + abs(((sm1(f1(:,1),3)+sm1(f1(:,2),3)+sm1(f1(:,3),3)))/3 - 1));
p1 = f1(NorthTri1,1);
p2 = f1(NorthTri1,2);
p3 = f1(NorthTri1,3);
if ~inpolygon(0,0,sm1(f1(NorthTri1,:),1),sm1(f1(NorthTri1,:),2))
    i = 1;
    while true
        if inpolygon(0,0,sm1(f1(i,:),1),sm1(f1(i,:),2)) && sum(sm1(f1(i,:),3)>0)==3
            NorthTri1 = i;
            p1 = f1(NorthTri1,1) ;p2 = f1(NorthTri1,2); p3 = f1(NorthTri1,3);
            break;
        else
            i = i+1;
        end
    end
end
[~,SouthTri1] = min(abs(((sm1(f1(:,1),1)+sm1(f1(:,2),1)+sm1(f1(:,3),1)))/3) + abs(((sm1(f1(:,1),2)+sm1(f1(:,2),2)+sm1(f1(:,3),2)))/3) + abs(((sm1(f1(:,1),3)+sm1(f1(:,2),3)+sm1(f1(:,3),3)))/3 + 1));
x1 = sm1(:,1)./(1-sm1(:,3)); 
y1 = sm1(:,2)./(1-sm1(:,3));
x2 = s2(:,1)./(1-s2(:,3)); 
y2 = s2(:,2)./(1-s2(:,3));
triedge_midptx = [mean(x1([p2 p3])), mean(x1([p1 p3])),mean(x1([p1 p2]))];
triedge_midpty = [mean(y1([p2 p3])), mean(y1([p1 p3])),mean(y1([p1 p2]))];
int_meanx = mean(x1(setxor(1:length(sm1), [p1,p2,p3])));
int_meany = mean(y1(setxor(1:length(sm1), [p1,p2,p3])));
[~, i] = max((triedge_midptx-int_meanx).^2 + (triedge_midpty-int_meany).^2);
x1(f1(NorthTri1,i)) = 3*int_meanx - sum(x1(f1(NorthTri1,setxor(1:3,i))));
y1(f1(NorthTri1,i)) = 3*int_meany - sum(y1(f1(NorthTri1,setxor(1:3,i))));

L = cotangent_laplacian([x1,y1],f1)/4;
M = L - sparse(lm1a, lm1a, lambda*ones(1, length(lm1a)), numofvertex1, numofvertex1);
c = zeros(numofvertex1,1); 
d = c;
c(lm1a) = -lambda*(x2(lm2a));
d(lm1a) = -lambda*(y2(lm2a));

fixed = [p1,p2,p3,find(~inpolygon(x1,y1,x1([p1 p2 p3]),y1([p1 p2 p3])))'];
[mrow,mcol,mval] = find(M(fixed,:));
M = M - sparse(fixed(mrow),mcol,mval,length(x1), length(x1)) + ...
        sparse(fixed,fixed,ones(size(fixed)),length(x1), length(x1));
c(fixed) = x1(fixed);
d(fixed) = y1(fixed);

z = M\complex(c,d);
x = real(z);
y = imag(z);
S = [2.*x./(1+x.^2+y.^2),2.*y./(1+x.^2+y.^2),(-1+x.^2+y.^2)./(1+x.^2+y.^2)];
fprintf('Done.\n');

%% Fix overlap
fprintf('Fixing overlap... \n');
x = S(:,1)./(1-S(:,3));
y = S(:,2)./(1-S(:,3));
x2 = s2(:,1)./(1-s2(:,3)); 
y2 = s2(:,2)./(1-s2(:,3));
v = [x1,y1];
f = f1([1:bigtri1-1,bigtri1+1:end],:);
initial_map = [x,y];
lm_id = lm1a';
lm_target = [x2(lm2a'),y2(lm2a')];
bdy_id = fixed';
bdy_target = initial_map(bdy_id,:);
Operator = createOperator(v',f');
map = initial_map;
alpha = 1;
beta = 1;
penalty = 0.01;
sigmaIncrease = 0.5;
mu_bound = 0.99;
balancing_ratio = exp(-sqrt(lambda)/4); % 1 if lambda = 0, 0 if lambda = inf

update_mu = beltrami_coefficient(v, f, initial_map);
overlap_count = sum(abs(update_mu)>=1);
mu_diff = max(abs(update_mu - zeros(length(f),1)));
landmark_error = mean(sqrt(sum((initial_map(lm_id,1:2)-lm_target).^2,2)));
fprintf('  Iteration  Mu_Difference  Overlap   Landmark error\n');
IterationNumber = 0;
    
fprintf('%7.0f %15f %8.0f %17f \n',[IterationNumber,mu_diff,overlap_count,landmark_error]);
    
while overlap_count ~= 0
    
    mu = update_mu;
    
    % Smooth BC
    penalty = penalty + sigmaIncrease;
    Smooth_Operator = speye(length(v)) + ...
        1/penalty*(alpha*speye(length(v)) - beta*L/2);
    smooth_mu = smoothing(update_mu,Smooth_Operator,Operator);
    smooth_mu(abs(smooth_mu)>=mu_bound) = ...
        smooth_mu(abs(smooth_mu)>=mu_bound)./...
        abs(smooth_mu(abs(smooth_mu)>=mu_bound))*mu_bound;

    % find BC direction for non overlap exact landmark matching
    map = linear_beltrami_solver(v,f,smooth_mu,[bdy_id;lm_id],[bdy_target;lm_target]); 
    exact_mu = beltrami_coefficient(v, f, map);
    
    % combine both direction, update BC
    combined_mu = exact_mu + balancing_ratio*(smooth_mu - exact_mu);
    combined_mu(abs(combined_mu)>=mu_bound) = ...
        combined_mu(abs(combined_mu)>=mu_bound)./...
        abs(combined_mu(abs(combined_mu)>=mu_bound))*mu_bound;
    
    % reconstruct without any landmark constraints
    map = linear_beltrami_solver(v,f,combined_mu,bdy_id,bdy_target); 
    update_mu = beltrami_coefficient(v, f, map);
        
    % evaluate the result
    landmark_error = mean(sqrt(sum((map(lm_id,1:2)-lm_target).^2,2)));
    mu_diff = max(abs(update_mu - mu));
    overlap_count = sum(abs(update_mu)>=1);
    IterationNumber = IterationNumber + 1;
    fprintf('%7.0f %15f %8.0f %17f \n',[IterationNumber,mu_diff,overlap_count,landmark_error]);
    
    if IterationNumber > 40
        warning('Iteration exceeds 40. Automatically terminated.');
        % consider changing the parameters to improve the performance
        break;
    end
end

x = map(:,1);
y = map(:,2);
S = [2.*x./(1+x.^2+y.^2),2.*y./(1+x.^2+y.^2),(-1+x.^2+y.^2)./(1+x.^2+y.^2)];
fprintf('Done.\n');

%% south pole LBS
DistfromSP = S(lm1a(unique(sort([1,index,index(length(index)-1)+1]))),3);
DistfromSP = sortrows([(1:length(DistfromSP))',DistfromSP],2); 
choice = 1; 
i = 2;
while(length(choice)~=2)
    if isempty(find(mod(DistfromSP(choice,1)-DistfromSP(i,1),length(lm1))==0, 1))
        choice = [choice,i];
    end
    i = i + 1;
end
new_south_centre = sum(S(f1(SouthTri1,1:3),1:3))/3;
new_south_centre = new_south_centre/norm(new_south_centre,2);
south_mobius_point = [S(lm1_endpt(DistfromNP(choice(1))),1:3);...
    S(lm1_endpt(DistfromNP(choice(2))),1:3);new_south_centre];
south_mobius_target = [s2(lm2_endpt(DistfromNP(choice(1))),1:3);...
    s2(lm2_endpt(DistfromNP(choice(2))),1:3);[0,0,-1]];
S = mobius_transformation(S,south_mobius_point,south_mobius_target);
fixed = find(S(:,3)<=max(S(lm1a,3)));
P = [S(:,1)./(1+S(:,3)), S(:,2)./(1+S(:,3))];
mu = beltrami_coefficient(P, f1, v0);
map = linear_beltrami_solver(P,f1,mu,fixed,P(fixed,1:2));
z = complex(map(:,1),map(:,2));

%% Final result
S_after = [2*real(z)./(1+abs(z).^2), 2*imag(z)./(1+abs(z).^2), -(abs(z).^2-1)./(1+abs(z).^2)];
map = mobius_transformation(S_after,south_mobius_target,south_mobius_point);
fprintf('Final landmark mismatch energy = %d\n', ...
    sum((s2(lm2a,1)-map(lm1a,1)).^2+(s2(lm2a,2)-map(lm1a,2)).^2 + (s2(lm2a,3)-map(lm1a,3)).^2));
