%% interpolate from 1 zgrid to another

grdFile0='/home/aruba/vol2/ipasmans/Exp/Prm/grd_ow2km_r13lin_mix.nc'; 
grdFile1='/home/aruba/vol2/ipasmans/Exp/Prm/grd_ow2km_r13lin_mix.nc';
dataFileIn='/home/aruba/vol2/ipasmans/Exp/Exp26f_2011_free/ocean_his_2388.nc';
dataFileOut='/home/aruba/vol2/ipasmans/Exp/Exp26f_2011_free/ini20110718.nc'; 
itIn=24;
itOut=1; 

%% read grid

grd0.h=ncread(grdFile0,'h'); 
grd0.maskr=ncread(grdFile0,'mask_rho'); 
grd0.maskv=ncread(grdFile0,'mask_v'); 
grd0.masku=ncread(grdFile0,'mask_u'); 

grd1.h=ncread(grdFile1,'h'); 

%% z

zeta=ncread(dataFileIn,'zeta',[1,1,itIn],[Inf,Inf,1]); 
sigma=[-39.5:-.5]/40;
z0=roms_sigma2z(sigma,grd0.h,zeta); 
z1=roms_sigma2z(sigma,grd1.h,zeta); 

%% interpolate

%time
val=ncread(dataFileIn,'ocean_time',[itIn],[1]); 
ncwrite(dataFileOut,'ocean_time',val,[itOut]); 
display( datestr(val/24/3600+datenum('2005-01-01')) ); 

%temp
display('temp')
val=ncread(dataFileIn,'temp',[1,1,1,itIn],[Inf,Inf,Inf,1]);
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
ncwrite(dataFileOut,'temp',val,[1,1,1,itOut]); 

%salt
display('salt')
val=ncread(dataFileIn,'salt',[1,1,1,itIn],[Inf,Inf,Inf,1]);
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
ncwrite(dataFileOut,'salt',val,[1,1,1,itOut]); 

%u
display('u')
z0u=.5*z0(1:end-1,:,:)+.5*z0(2:end,:,:); 
z1u=.5*z1(1:end-1,:,:)+.5*z1(2:end,:,:); 
val=ncread(dataFileIn,'u',[1,1,1,itIn],[Inf,Inf,Inf,1]);
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
ncwrite(dataFileOut,'u',val,[1,1,1,itOut]); 

%v
display('v')
z0v=.5*z0(:,1:end-1,:)+.5*z0(:,2:end,:); 
z1v=.5*z1(:,1:end-1,:)+.5*z1(:,2:end,:); 
val=ncread(dataFileIn,'v',[1,1,1,itIn],[Inf,Inf,Inf,1]);
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
ncwrite(dataFileOut,'v',val,[1,1,1,itOut]); 

%zeta
display('zeta')
val=ncread(dataFileIn,'zeta',[1,1,itIn],[Inf,Inf,1]); 
val(isnan(val))=0; 
ncwrite(dataFileOut,'zeta',val,[1,1,itOut]); 

%ubar
display('ubar/vbar')
val=ncread(dataFileIn,'ubar',[1,1,itIn],[Inf,Inf,1]); 
val(isnan(val))=0; 
ncwrite(dataFileOut,'ubar',val,[1,1,itOut]); 

val=ncread(dataFileIn,'vbar',[1,1,itIn],[Inf,Inf,1]); 
val(isnan(val))=0; 
ncwrite(dataFileOut,'vbar',val,[1,1,itOut]); 



