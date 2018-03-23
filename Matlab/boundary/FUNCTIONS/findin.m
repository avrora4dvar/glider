function [IN]=findin(V,LIMS);

% AK, 6/21/04
% FINDIN.M finds indices (IN) of elements of a 1D vector V in the interval LIMS
% USAGE:
%
% [IN]=findin(V,LIMS);
%
% LIMS=[x1 x2];


nv=size(V);             % <- check V is 1D
if length(nv)>2 
 error('in findin: dim(V)>2');
elseif prod(nv)>max(nv)
 error('in findin: V is not 1D');
end

IN=find(V>=LIMS(1) & V<=LIMS(2));

