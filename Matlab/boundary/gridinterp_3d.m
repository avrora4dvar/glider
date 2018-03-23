function [a]=gridinterp_3d(XX,YY,z,A,in0,iclosest,x,y,zroms)

% extend values to lower NaN hycom layers
[nx,ny,nz]=size(A);
A1=A;
mask3d=ones([nx ny nz]);
mask3d(isnan(A))=0;
[tmp,ic]=max(mask3d,[],3); % ic is a 2d array of lowest ~NaN hycom cells
for l=1:ny
 for k=1:nx
  if mask3d(k,l,end)==1
   A1(k,l,1:ic(k,l)-1)=A1(k,l,ic(k,l));
  end
 end
end

% fill in masked values (use info on mask and closest unmasked neighbor
% from SSH)
if ~isempty(in0)
 for k=1:nz
  tmp=A1(:,:,k);
  tmp(in0)=tmp(iclosest);
  A1(:,:,k)=tmp;
 end
end

Y3D=repmat(YY,[1 1 nz]);
X3D=repmat(XX,[1 1 nz]);
Z3D=repmat(reshape(z,[1 1 nz]),[nx ny 1]); 

a=interp3(Y3D,X3D,Z3D,A1,y,x,zroms);