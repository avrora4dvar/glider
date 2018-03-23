function [ out ] = get_roms_cross( varargin )
%GET_ROMS_CROSSget model results along zonal cross-section
%   get_roms_cross_zonal(grid_file,roms_dir,field_name,grd,loc)
%   gets the zonal cross-section along lat=loc.lat between min(loc.lon) and
%   max(loc.lon) at times loc.t) for the variable field_name defined on the
%   grid in grd. Here grid_file is the grid file for the model and roms_dir
%   is the directory containing the history files. 

%input
cross_type=varargin{1};
grid_file=varargin{2}; 
roms_dir=varargin{3}; 
field_name=varargin{4}; 
grd=varargin{5};
loc=varargin{6};
if length(varargin)==7
    filter=varargin{7};
else
    filter=false; 
end

%s-grid parameters
par.Vtransform=2; 
par.Vstretching=4;
par.Tcline=50; %[m]
par.theta_s=8;
par.theta_b=3;
par.n=40; 
par.hc=par.Tcline; %for Vtransform 1: par.hc=min(par.Tcline,min(grd0.h(:))); 

%grid indices
ilat(1)=find(grd.lat(1,:)<=min(loc.lat),1,'last'); 
ilat(2)=find(grd.lat(1,:)>=max(loc.lat),1,'first'); 
if ilat(1)==ilat(2)
    wlat=.5; 
else
    wlat=(loc.lat-grd.lat(1,ilat(1)))/(grd.lat(1,ilat(2))-grd.lat(1,ilat(1))); 
end
ilon(1)=find(grd.lon(:,1)<=min(loc.lon),1,'last'); 
ilon(2)=find(grd.lon(:,1)>=max(loc.lon),1,'first'); 
if ilon(1)==ilon(2)
    wlon=.5; 
else
    wlon=(loc.lon-grd.lon(ilon(1),1))/(grd.lon(ilon(2),1)-grd.lon(ilon(1),1)); 
end

%bottom along cross-section
h.lon=ncread(grid_file,'lon_rho'); 
h.lat=ncread(grid_file,'lat_rho'); 
if length(loc.lat)==1
    h.val=(1-wlat)*ncread(grid_file,'h',[1,ilat(1)],[Inf,1])+...
        wlat*ncread(grid_file,'h',[1,ilat(2)],[Inf,1]);
    h.val=interp1(h.lon(:,1),h.val,grd.lon(ilon(1):ilon(2),1));
elseif length(loc.lon)==1
    h.val=(1-wlon)*ncread(grid_file,'h',[ilon(1),1],[1,Inf])+...
        wlon*ncread(grid_file,'h',[ilon(2),1],[1,Inf]);
    h.val=interp1(h.lat(1,:),h.val,grd.lat(1,ilat(1):ilat(2)));
else
    error('loc.lat or loc.lon must have length of 1');
end

%get roms indices
if filter
    tLim=[min(loc.t),max(loc.t)]+[-61,61]/24;
else
    tLim=[min(loc.t),max(loc.t)]+[-1,1]/24;
end
[hash,t_roms]=get_hash_roms(roms_dir,tLim);


%read roms\
if length(loc.lat)==1
    info=ncinfo(hash(1).filename,field_name);
    val1=nan(diff(ilon)+1,info.Size(3)+2,length(t_roms));
    zeta1=nan(diff(ilon)+1,length(t_roms));
    iBnd=[1,0];
    for i0=1:length(hash)
        iBnd(2)=iBnd(2)+length(hash(i0).index);
        %field
        val1(:,2:end-1,iBnd(1):iBnd(2))=squeeze(...
            (1-wlat)*ncread(hash(i0).filename,field_name,[ilon(1),ilat(1),1,hash(i0).index(1)],...
            [diff(ilon)+1,1,Inf,length(hash(i0).index)])+wlat*...
            ncread(hash(i0).filename,field_name,[ilon(1),ilat(2),1,hash(i0).index(1)],...
            [diff(ilon)+1,1,Inf,length(hash(i0).index)]) );
        %zeta
        zeta2=squeeze(...
            (1-wlat)*ncread(hash(i0).filename,'zeta',[1,ilat(1),hash(i0).index(1)],...
            [Inf,1,length(hash(i0).index)])+wlat*...
            ncread(hash(i0).filename,'zeta',[1,ilat(2),hash(i0).index(1)],...
            [Inf,1,length(hash(i0).index)]) );
        zeta1(:,iBnd(1):iBnd(2))=interp1(h.lon(:,1),zeta2,grd.lon(ilon(1):ilon(2),1));
        iBnd(1)=iBnd(2)+1;
    end
    val1(:,1,:)=val1(:,2,:); val1(:,end,:)=val1(:,end-1,:);
elseif length(loc.lon)==1
    info=ncinfo(hash(1).filename,field_name);
    val1=nan(diff(ilat)+1,info.Size(3)+2,length(t_roms));
    zeta1=nan(diff(ilat)+1,length(t_roms));
    iBnd=[1,0];
    for i0=1:length(hash)
        iBnd(2)=iBnd(2)+length(hash(i0).index);
        %field
        val1(:,2:end-1,iBnd(1):iBnd(2))=squeeze(...
            (1-wlon)*ncread(hash(i0).filename,field_name,[ilon(1),ilat(1),1,hash(i0).index(1)],...
            [1,diff(ilat)+1,Inf,length(hash(i0).index)])+wlon*...
            ncread(hash(i0).filename,field_name,[ilon(2),ilat(1),1,hash(i0).index(1)],...
            [1,diff(ilat)+1,Inf,length(hash(i0).index)]) );
        %zeta
        zeta2=squeeze(...
            (1-wlon)*ncread(hash(i0).filename,'zeta',[ilon(1),1,hash(i0).index(1)],...
            [1,Inf,length(hash(i0).index)])+wlon*...
            ncread(hash(i0).filename,'zeta',[ilon(2),1,hash(i0).index(1)],...
            [1,Inf,length(hash(i0).index)]) );
        zeta1(:,iBnd(1):iBnd(2))=interp1(h.lat(1,:),zeta2,grd.lat(1,ilat(1):ilat(2)));
        iBnd(1)=iBnd(2)+1;
    end
    val1(:,1,:)=val1(:,2,:); val1(:,end,:)=val1(:,end-1,:);
end

%filter
if filter
    val1=filter_lanczos(val1,1,1/40,60,3); 
    zeta=filter_lanczos(zeta,1,1/40,60,2); 
end

%interpolate in time
[tUni,t2uni,uni2t]=unique(loc.t); 
val1=permute(val1,[3,1,2]); 
zeta1=permute(zeta1,[2,1]); 
if length(t_roms)>1
    val1=interp1(t_roms,val1,tUni);
    zeta1=interp1(t_roms,zeta1,tUni); 
end
val1=permute(val1,[2,3,1]); 
zeta1=permute(zeta1,[2,1]); 
val=val1(:,:,uni2t); 
zeta=zeta1(:,uni2t); 

%calculate z
h.val=repmat(h.val,[1,size(zeta,2)]); 
s=[-par.n,-par.n+.5:1:-.5,0]/par.n; 
z=z24(s,zeta(:),h.val(:),par);
z=reshape(z,[size(z,1),size(val,1),size(val,3)]); 
z=permute(z,[2,1,3]); 

if length(loc.lat)==1
    out=struct('lon',repmat(grd.lon(ilon(1):ilon(2),1),[1,length(s)]),'lat',loc.lat(1),...
        'z',z,'t',loc.t,'field',field_name,'val',val);
else
    out=struct('lat',repmat(grd.lat(1,ilat(1):ilat(2))',[1,length(s)]),'lon',loc.lon(1),...
        'z',z,'t',loc.t,'field',field_name,'val',val);
end

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




