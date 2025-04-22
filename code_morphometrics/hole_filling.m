function f_filled = hole_filling(f,n)

% Fill unwanted holes in a triangulation.
%
% f: nf x 3 triangulation
% n: number of holes to keep 
%    (for genus-0 closed surfaces, n should be 0)
%    (for simply-connected open surfaces, n should be 1)
%
% Written by Gary Choi, 2023

bd = meshboundaries(f);

if length(bd) < n
    disp('Actual number of holes < desired number of holes!');
    f_filled = f;
    return;
end

% keep the first n holes unfilled, fill the remaining holes
f_filled = f;
for i = (n+1):length(bd)
    bdi = bd{i};
    k = floor(length(bdi)/2);
    
    f_temp = fliplr(bdi([(1:(k-1))', (2:k)', (length(bdi)+1-(1:(k-1)))'; ...
        (2:k)', (length(bdi)+1-(2:k))', (length(bdi)+1-(1:(k-1)))']));
    
    f_filled = [f_filled; f_temp];
    
end
