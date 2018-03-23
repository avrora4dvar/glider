%calculate balance-operator covariance between one point and the rest
clear all; clc; 
addpath('../'); 

grdFile='/home/aruba/vol2/ipasmans/Exp/Prm/grd_ow2km_r13lin_mix.nc'; 
sigFile='/home/aruba/vol2/ipasmans/Exp/Prm/per_sig_bal.nc'; 
backFile='/home/aruba/vol2/ipasmans/Exp/Exp35/Exp35_2392/NL/ocean_his_0001.nc';
tmpFile='/home/aruba/vol1/ipasmans/Exp/tmp_cov'; 
binDir='/home/aruba/vol2/ipasmans/Bin_Exp35/';
outFile='/home/server/student/homes/ipasmans/Data/article_glider/cov_bal_v.mat'; 
figFile='/home/server/student/homes/ipasmans/Figures/article_glider/cov_bal_v'; 


loc=struct('lon', -124.3,'lat',44.642,'field','temp','z',0); 
loc=struct('lon', -124.539,'lat',44.684,'field','temp','z',0);  %mean glider transect
loc=struct('lon', -127.539,'lat',44.684,'field','v','z',0);  %mean glider transect
R0=1025; Tcoef=1.7e-4; Scoef=7.5e-4; 
dl=[25e3,25e3,15]; 
sig=.9; z0=100;  

%exp35: cov referenced in point 203,132 h=250m


%% read grid


nLayer=ncinfo(backFile,'temp');
nLayer=nLayer.Size(3); 
[grdr,grdu,grdv]=read_grid_file(grdFile,nLayer); 




[~,loc.ix]=min(abs(grdr.lon(:,1)-loc.lon)); 
[~,loc.iy]=min(abs(grdr.lat(1,:)-loc.lat));
[~,loc.iz]=min(abs(grdr.zr(loc.ix,loc.iy,:)-loc.z)); 

if strcmpi(loc.field,'u')
    grd=grdu; 
elseif strcmpi(loc.field,'v')
    grd=grdv; 
else
    grd=grdr; 
end

[~,grd.ilat]=min(abs(grd.lat(1,:)-loc.lat)); 


ix(1)=find(grd.lon(:,1)<=loc.lon,1,'last'); 
ix(2)=find(grd.lon(:,1)>loc.lon,1,'first');
wx(1)=grd.lon(ix(2),1)-loc.lon;
wx(2)=loc.lon-grd.lon(ix(1),1); 
wx=wx/sum(wx); 

iy(1)=find(grd.lat(1,:)<=loc.lat,1,'last'); 
iy(2)=find(grd.lat(1,:)>loc.lat,1,'first'); 
wy(1)=grd.lat(1,iy(2))-loc.lat;
wy(2)=loc.lat-grd.lat(1,iy(1)); 
wy=wy/sum(wy); 

if ~strcmpi(loc.field,'zeta')
    z1=zeros(1,1,nLayer);
    for i1=1:2
        for i2=1:2
            z1=z1+grd.zr(ix(i1),iy(i2),:)*wx(i1)*wy(i2);
        end
    end
    z1=squeeze(z1);
    if loc.z>=max(z1)
        iz(1)=nLayer; iz(2)=nLayer;
        wz=[.5,.5];
    elseif loc.z<=min(z1)
        iz(1)=1; iz(2)=1;
        wz=[.5,.5];
    else
        iz(1)=find(z1<=loc.z,1,'last');
        iz(2)=find(z1>loc.z,1,'first');
        wz(1)=z1(iz(2))-loc.z;
        wz(2)=loc.z-z1(iz(1));
        wz=wz/sum(wz);
    end
end


if true
    
    copyfile(sigFile,[tmpFile,'_in.nc']); 
%% create input

val=ncread([tmpFile,'_in.nc'],'temp'); val=zeros(size(val)); 
if strcmpi(loc.field,'temp')
    for i1=1:2
        for i2=1:2
            for i3=1:2
                val(ix(i1),iy(i2),iz(i3))=val(ix(i1),iy(i2),iz(i3))+...
                    wx(i1)*wy(i2)*wz(i3); 
            end
        end
    end
end
ncwrite([tmpFile,'_in.nc'],'temp',val); 

val=ncread([tmpFile,'_in.nc'],'salt'); val=zeros(size(val)); 
if strcmpi(loc.field,'salt')
    for i1=1:2
        for i2=1:2
            for i3=1:2
                val(ix(i1),iy(i2),iz(i3))=val(ix(i1),iy(i2),iz(i3))+...
                    wx(i1)*wy(i2)*wz(i3);
            end
        end
    end
end
ncwrite([tmpFile,'_in.nc'],'salt',val); 

val=ncread([tmpFile,'_in.nc'],'u'); val=zeros(size(val)); 
if strcmpi(loc.field,'u')
    for i1=1:2
        for i2=1:2
            for i3=1:2
                val(ix(i1),iy(i2),iz(i3))=val(ix(i1),iy(i2),iz(i3))+...
                    wx(i1)*wy(i2)*wz(i3);
            end
        end
    end 
end
ncwrite([tmpFile,'_in.nc'],'u',val); 

val=ncread([tmpFile,'_in.nc'],'v'); val=zeros(size(val)); 
if strcmpi(loc.field,'v')
    for i1=1:2
        for i2=1:2
            for i3=1:2
                val(ix(i1),iy(i2),iz(i3))=val(ix(i1),iy(i2),iz(i3))+...
                    wx(i1)*wy(i2)*wz(i3);
            end
        end
    end
end
ncwrite([tmpFile,'_in.nc'],'v',val); 

val=ncread([tmpFile,'_in.nc'],'zeta'); val=zeros(size(val)); 
if strcmpi(loc.field,'zeta')
    for i1=1:2
        for i2=1:2
            val(ix(i1),iy(i2))=val(ix(i1),iy(i2))+...
                    wx(i1)*wy(i2);
        end
    end
end
ncwrite([tmpFile,'_in.nc'],'zeta',val); 


%% ad_balance_3D

display('ad_balance'); 
fid=fopen([tmpFile,'.in'],'w+'); 
fprintf(fid,'#\n%s\n',grdFile); 
fprintf(fid,'#\n%s\n',[tmpFile,'1.nc']); 
fprintf(fid,'#\n%s\n',[tmpFile,'_in.nc']); 
fprintf(fid,'#\n%s\n',backFile); 
fprintf(fid,'#\n%f,%f,%f\n',R0,Tcoef,Scoef);
fclose(fid); 

command=sprintf('%s/ad_balance_3D < %s',binDir,[tmpFile,'.in']); 
unix(command); 



%% cov_balance

display('cov_balance'); 
fid=fopen([tmpFile,'.in'],'w+'); 
fprintf(fid,'#\n%s\n',grdFile); 
fprintf(fid,'#\n%s\n',[tmpFile,'1.nc']); 
fprintf(fid,'#\n%f\n',dl(1)); 
fprintf(fid,'#\n%f\n',dl(2)); 
fprintf(fid,'#\n%f\n',dl(3)); 
fprintf(fid,'#\n%s\n',sigFile); 
fprintf(fid,'#\n%f\n',z0); 
fid=fclose(fid); 

command=sprintf('%s/cov_bal_1_1 < %s',binDir,[tmpFile,'.in']); 
unix(command); 



%% tl_balance_3D

display('tl_balance'); 
fid=fopen([tmpFile,'.in'],'w+'); 
fprintf(fid,'#\n%s\n',grdFile); 
fprintf(fid,'#\n%s\n',[tmpFile,'1.nc']); 
fprintf(fid,'#\n%s\n',[tmpFile,'_out.nc']); 
fprintf(fid,'#\n%s\n',''); 
fprintf(fid,'#\n%f,%f,%f\n',R0,Tcoef,Scoef);
fclose(fid); 

command=sprintf('%s/tl_balance_3D < %s',binDir,[tmpFile,'.in']); 
unix(command); 

%% extract and save

cov.temp=ncread([tmpFile,'_out.nc'],'temp'); 
cov.salt=ncread([tmpFile,'_out.nc'],'salt'); 
cov.u=ncread([tmpFile,'_out.nc'],'u'); 
cov.v=ncread([tmpFile,'_out.nc'],'v'); 
cov.zeta=ncread([tmpFile,'_out.nc'],'zeta'); 

delete([tmpFile,'1.nc']); 
delete([tmpFile,'_in.nc']); 
delete([tmpFile,'_out.nc']); 
delete([tmpFile,'.in']); 

save(outFile,'cov','loc','grd'); 
end

%% Print

flag_plot=false; 
if flag_plot
   load(outFile); 
   
   cov.temp=squeeze(cov.temp(:,:,end,1)); 
   cov.salt=squeeze(cov.salt(:,:,end,1));
   cov.u=squeeze(cov.u(:,:,end,1)); 
   cov.v=squeeze(cov.v(:,:,end,1)); 
   
   
    
   figure(); set(gcf,'color','w','units','inches','renderer','zbuffer'); hold on; 
   set(gcf,'papersize',[5.5,5.5],'position',[0,0,5.5,5.5]); 
   
   
   ctickT=[-.06:.01:.06]; ctickS=[-.15:.01:.15]; 
   contourf(grdr.lon,grdr.lat,cov.salt,ctickS,'linestyle','none'); 
   
   contour(grdr.lon,grdr.lat,cov.temp,ctickT(ctickT>0),'k','linewidth',1.5,'linestyle','-');  
   contour(grdr.lon,grdr.lat,cov.temp,ctickT(ctickT<0),'k','linewidth',1.5,'linestyle','--');  
   
   step=3; s=10; 
   u1=.5*cov.u(1:end-1,2:end-1)+.5*cov.u(2:end,2:end-1); 
   v1=.5*cov.v(2:end-1,1:end-1)+.5*cov.v(2:end-1,2:end); 
   in=hypot(u1,v1)<.1*max(hypot(u1(:),v1(:))); u1(in)=NaN; v1(in)=NaN; 
   quiver(grd.lon(2:step:end-1,2:step:end-1),grd.lat(2:step:end-1,2:step:end-1),...
      s*u1(1:step:end,1:step:end)./cos(deg2rad(grd.lat(2:step:end-1,2:step:end-1))),...
      s*v1(1:step:end,1:step:end),0,'color',[1,1,1]*.1,'linewidth',1); 
  
   contour(grd.lon,grd.lat,grd.h,[200,1000],'linestyle',':','color',[1,1,1]*0,'linewidth',1.5); 
   h=plot(loc.lon,loc.lat,'kx','markersize',8); 
   
   p=plot(loc.lon,loc.lat,'mx','markersize',12,'linewidth',2);
   
   %add coastline
   load('../coastLine.mat'); 
   set(gca,'layer','top');
   for k=1:length(coastLine)
       patch(coastLine(k).Lon(1:end-1),coastLine(k).Lat(1:end-1),...
           ones(size(coastLine(k).Lat(1:end-1)))*1,[.5,.6,.5]);
   end
   
   cmap=jet(length(ctickS)-1); 
   cmap(round(.5*size(cmap,1)+[0:1]),:)=repmat([1,1,1],[2,1]); 
   colormap(cmap)
   set(gca,'clim',[min(ctickS),max(ctickS)],'xlim',[-1,1]+loc.lon,...
       'ylim',[-1,1]+loc.lat,'dataaspectratio',[1,cos(deg2rad(loc.lat)),1]);
   set(gca,'layer','top','box','on'); 
   hCB=colorbar(); set(get(hCB,'title'),'string','[C ppt]'); 
   
   print(gcf,sprintf('%s_%s_lay%-.2d',figFile,loc.field,loc.iz),'-r200','-dpng'); 
   
end


