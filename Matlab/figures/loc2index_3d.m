function [ix,iy,iz,w,loc] = loc2index_3d(grid_file,roms_dir,grd,loc)
%LOC2INDEX_Z [ix,iy,iz,w,loc]=loc2index_3d(grid_file,roms_dir,grd,loc)
%Calculates indices and weights or 3D grid

npoints=length(loc.lon);
nt=length(loc.t); 

%rho grid
grd0.lon=ncread(grid_file,'lon_rho'); 
grd0.lat=ncread(grid_file,'lat_rho'); 
grd0.h=ncread(grid_file,'h'); 

%s-grid parameters
par.Vtransform=2; 
par.Vstretching=4;
par.Tcline=50; %[m]
par.theta_s=8;
par.theta_b=3;
par.n=40; 
par.hc=par.Tcline; %for Vtransform 1: par.hc=min(par.Tcline,min(grd0.h(:))); 

%interpolate h
h=interp2(grd0.lon',grd0.lat',grd0.h',loc.lon,loc.lat); 

%get zeta
zeta=get_roms_2d(roms_dir,'zeta',grd,loc);
zeta=zeta.val; 


%perform depth/z conversion
if isfield(loc,'z') && isfield(loc,'depth')
    error('loc must contain either field z or depth, not both'); 
elseif isfield(loc,'z')
    loc.depth=zeta-loc.z;
elseif isfield(loc,'depth')
    loc.z=zeta-loc.depth;
else
    error('loc must at least contain either field z or depth'); 
end

%construct z-grid
s=[-par.n+.5:-.5]/par.n; 
z=z24(s,zeta,h,par); 

%horizontal index
[ix,iy,w]=loc2index_2d(grd,loc); 
ix=repmat(ix,[2,1]); iy=repmat(iy,[2,1]); w=repmat(w,[2,1]); 
iz=nan(size(ix)); 

%vertical index
for i2=1:npoints    
   if any(isnan([loc.lon(i2),loc.lat(i2),loc.z(i2)]))
       w(:,i2)=0; ix(:,i2)=NaN; iy(:,i2)=NaN; iz(:,i2)=NaN; 
       continue;
   end
   if all(w(:,i2)==0); continue; end
   if loc.z(i2)<=z(1,i2)
       iz(5:8,i2)=1; 
       w(1:4,i2)=0; 
   elseif loc.z(i2)>=z(end,i2)
       iz(1:4,i2)=size(z,1); 
       w(5:8,i2)=0; 
   else
       iz1=interp1(reshape(z(:,i2),1,[]),[1:size(z,1)],loc.z(i2)); 
       iz(1:4,i2)=floor(iz1); 
       w(1:4,i2)=w(1:4,i2)*(1-abs(iz1-iz(1,i2))); 
       iz(5:8,i2)=ceil(iz1); 
       w(5:8,i2)=w(5:8,i2)*(1-abs(iz1-iz(5,i2))); 
   end     
end

%Normalize w
sumw=sum(w,1); 
in=sumw>0;
w(:,in)=bsxfun(@times,w(:,in),1./sumw(in)); 

end

%%


function out=C4(s,par)
%ROMS C-function for Vstretching=4

out=(1-cosh(par.theta_s*s))/(cosh(par.theta_s)-1); 
out=(exp(par.theta_b*out)-1)/(1-exp(-par.theta_b)); 

end

%% 

function out=S2(s,h,par)
%ROMS S-function for Vtransform=2

h=h(:)'; 
s=s(:);

C=C4(s,par); 

h=repmat(h,[size(s,1),1]); 
C=repmat(C,[1,size(h,2)]);
s=repmat(s,[1,size(h,2)]); 

out=(par.hc*s+h.*C)./(par.hc+h); 

end

%% 

function z=z24(s,zeta,h,par)
%ROMS s-z conversion function for Vtransform=2,Vstretching=4

s=reshape(s,[],1); 
h=reshape(h,1,[]); 
zeta=reshape(zeta,1,[]); 

s=repmat(s,[1,size(h,2)]); 
zeta=repmat(zeta,[size(s,1),1]); 
h=repmat(h,[size(s,1),1]); 

z=zeta+(zeta+h).*S2(s(:,1),h(1,:),par); 

end




