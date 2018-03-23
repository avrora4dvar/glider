function [speed,deg] = cartesian2nautical(u,v,rot,dim)
%CARTESIAN2NAUTICAL 
% [speed,deg]=cartesian2nautical(u,v,rot,dim) convert the velocity given with respect to
% the east,north vector to speeds and degrees North. dim is the index of
% the time dimension. rot is fixed rotation

%speed
speed=hypot(u,v); 

%reshape matrix
dimList=[1:ndims(u)]; 
dimList=circshift(dimList,[0,-dim+1]); 
u=permute(u,dimList); 
v=permute(v,dimList); 
s=size(u); 
u=reshape(u,size(u,1),[]); 
v=reshape(v,size(v,1),[]); 

deg=nan(size(u)); 
for i2=1:size(u,2)
   uv=u(:,i2)+1i*v(:,i2); 
   inan=~isnan(uv) & abs(uv)~=0; 
   uv=uv(inan); 
   
   %rotate plane
   rot=mod(rot,360); 
   uv=uv*exp(1i*deg2rad( (rot-270) )); 
   theta=rad2deg(angle(uv));
   
   for i1=2:length(theta)
       dt=theta(i1)-theta(i1-1); 
       if dt>180; theta(i1)=theta(i1)-360; end
       if dt<180; theta(i1)=theta(i1)+360; end
   end
   
   theta=theta-(rot-270); 
   
   %cartesian to nautical
   theta=-theta+90;
   
   %mean between -180,180
   theta=theta-round(mean(theta)/360)*360;
   deg(inan,i2)=theta; 
end


%transform back
deg=reshape(deg,s); 
deg=ipermute(deg,dimList); 

end

