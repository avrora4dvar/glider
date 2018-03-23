%% interpolate from 1 zgrid to another

grdFile0='/home/aruba/vol2/ipasmans/Exp/Prm/grd_ow2km_r13lin_mix.nc'; 
grdFile1='/home/aruba/vol2/ipasmans/Exp/Prm/grd_ow4km_interp.nc';
dataFile='/home/aruba/vol2/ipasmans/Exp/Exp27_2344/Iter0/ini.nc';
it=1; 

%% read grid

grd0.h=ncread(grdFile0,'h'); 
grd0.lon=ncread(grdFile0,'lon_rho'); 
grd0.lat=ncread(grdFile0,'lat_rho'); 
grd0.maskr=ncread(grdFile0,'mask_rho'); 
grd0.maskv=ncread(grdFile0,'mask_v'); 
grd0.masku=ncread(grdFile0,'mask_u'); 

grd1.h=ncread(grdFile1,'h'); 
grd1.lon=ncread(grdFile1,'lon_rho'); 
grd1.lat=ncread(grdFile1,'lat_rho'); 



%% z

zeta=ncread(dataFile,'zeta',[1,1,it],[Inf,Inf,1]); 
sigma=[-39.5:-.5]/40;
z0=roms_sigma2z(sigma,grd0.h,zeta); 
z1=roms_sigma2z(sigma,grd1.h,zeta); 

%% interpolate

%temp
display('temp')
val=ncread(dataFile,'temp',[1,1,1,it],[Inf,Inf,Inf,1]);
for i1=1:size(grd0.maskr,1)
    for i2=1:size(grd0.maskr,2)
        if ~grd0.maskr(i1,i2); continue; end
        
        valloc=squeeze(val(i1,i2,:));
        z0loc=squeeze(z0(i1,i2,:));
        z1loc=squeeze(z1(i1,i2,:));
        
        valout=interp1(z0loc,valloc,z1loc);
        valout(z1loc>max(z0loc))=valloc(end);
        valout(z1loc<min(z0loc))=valloc(1);
        
        val(i1,i2,:)=valout;
        
    end
end
val(isnan(val))=10; 
ncwrite(dataFile,'temp',val,[1,1,1,it]); 

%salt
display('salt')
val=ncread(dataFile,'salt',[1,1,1,it],[Inf,Inf,Inf,1]);
for i1=1:size(grd0.maskr,1)
    for i2=1:size(grd0.maskr,2)
        if ~grd0.maskr(i1,i2); continue; end
        
        valloc=squeeze(val(i1,i2,:));
        z0loc=squeeze(z0(i1,i2,:));
        z1loc=squeeze(z1(i1,i2,:));
        
        valout=interp1(z0loc,valloc,z1loc);
        valout(z1loc>max(z0loc))=valloc(end);
        valout(z1loc<min(z0loc))=valloc(1);
        
        val(i1,i2,:)=valout;
        
    end
end
val(isnan(val))=34; 
ncwrite(dataFile,'salt',val,[1,1,1,it]); 

%u
display('u')
z0u=.5*z0(1:end-1,:,:)+.5*z0(2:end,:,:); 
z1u=.5*z1(1:end-1,:,:)+.5*z1(2:end,:,:); 
val=ncread(dataFile,'u',[1,1,1,it],[Inf,Inf,Inf,1]);
for i1=1:size(grd0.masku,1)
    for i2=1:size(grd0.masku,2)
        if ~grd0.masku(i1,i2); continue; end
        
        valloc=squeeze(val(i1,i2,:));
        z0loc=squeeze(z0u(i1,i2,:));
        z1loc=squeeze(z1u(i1,i2,:));
        
        valout=interp1(z0loc,valloc,z1loc);
        valout(z1loc>max(z0loc))=valloc(end);
        valout(z1loc<min(z0loc))=valloc(1);
        
        val(i1,i2,:)=valout;
        
    end
end
val(isnan(val))=0; 
ncwrite(dataFile,'u',val,[1,1,1,it]); 

%v
display('v')
z0v=.5*z0(:,1:end-1,:)+.5*z0(:,2:end,:); 
z1v=.5*z1(:,1:end-1,:)+.5*z1(:,2:end,:); 
val=ncread(dataFile,'v',[1,1,1,it],[Inf,Inf,Inf,1]);
for i1=1:size(grd0.maskv,1)
    for i2=1:size(grd0.maskv,2)
        if ~grd0.maskv(i1,i2); continue; end
        
        valloc=squeeze(val(i1,i2,:));
        z0loc=squeeze(z0v(i1,i2,:));
        z1loc=squeeze(z1v(i1,i2,:));
        
        valout=interp1(z0loc,valloc,z1loc);
        valout(z1loc>max(z0loc))=valloc(end);
        valout(z1loc<min(z0loc))=valloc(1);
        
        val(i1,i2,:)=valout;
        
    end
end
val(isnan(val))=0; 
ncwrite(dataFile,'v',val,[1,1,1,it]); 

%zeta
display('zeta')
val=ncread(dataFile,'zeta'); 
val(isnan(val))=0; 
ncwrite(dataFile,'zeta',val); 

%ubar
display('ubar/vbar')
val=ncread(dataFile,'ubar'); 
val(isnan(val))=0; 
ncwrite(dataFile,'ubar',val); 

val=ncread(dataFile,'vbar'); 
val(isnan(val))=0; 
ncwrite(dataFile,'vbar',val); 



