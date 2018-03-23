function [isoZ] = mesh_isosurface(z,val,isoVal)
%ISOPLANE Get z-coordinates of an isoplane
%   [isoZ]=isoplane(z,val,isoVal)
%   z: 3D-array with z-coordinates
%   val: 3D-array with values variable on z-coordinates
%   isoVal: 1D-array with values variable on isoplanes
%   isoZ: 3D-array with isoZ(:,:,k) the z-coordinates of the isoplane with 
%         value isoVal(k)

%initiate output
isoZ=nan(size(z,1),size(z,2),numel(isoVal)); 

for k=1:length(isoVal)
    for i2=1:size(z,2)
        for i1=1:size(z,1)
            
            if isempty(iL
                iL=find(val(i1,i2,:)<=isoVal(k),1,
            
            
        end
    end
end

end

