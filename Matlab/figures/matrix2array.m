function [array] = matrix2array(matrix,dim,s)
%MATRIX2ARRAY Transforms 2D-matrix to array
%   [array]=matrix2array(matrix,dim,s) transforms a 2D-matrix
%   to an array of size s. The first dimension in matrix becomes the dimth
%   dimension in the array. This function is the inverse of array2matrix

dimList=[1:length(s)]; 
dimList=circshift(dimList,[0,-dim+1]); 
s=circshift(s,[0,-dim+1]); 

matrix=reshape(matrix,s); 
array=ipermute(matrix,dimList); 

end

