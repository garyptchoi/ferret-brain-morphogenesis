function [vertex,face] = read_off(name)

if ~strcmp(name(end-3:end),'.off')
    name = [name,'.off'];
end
fid = fopen(name);
textscan(fid,'OFF%*[^\n]',1);
textscan(fid,'COFF%*[^\n]',1);
A = cell2mat(textscan(fid,'%f %f%*[^\n]',1));
numofvertex = A(1,1);
clear A;
vertex = cell2mat(textscan(fid,'%f%f%f%*[^\n]',numofvertex));
face = cell2mat(textscan(fid,'3 %f%f%f%*[^\n]'))+1;
fclose(fid);