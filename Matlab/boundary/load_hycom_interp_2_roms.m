function [SSH,U,V,Thycom,Shycom,ssh_roms,u_roms,v_roms,T_roms,S_roms]=...
     load_hycom_interp_2_roms...
     (url_2d,url_uvel,url_vvel,url_temp,url_salt,...
     XX,YY,z,i1,j1,in0,iclosest,it1,...
     xbr, ybr,xbr2,ybr2,xbu2,ybu2,xbv2,ybv2,...
     zr,zu,zv)
 
 [nx,ny]=size(XX);
 nz=length(z);
 
 % ssh:
 disp('ssh...');
 SSH=ncread(url_2d,'ssh',[i1 j1 it1],[nx ny 1]);
 ssh_roms=gridinterp_2d(XX,YY,SSH,in0,iclosest,xbr,ybr); 
 
 % u: x
 disp('u...');
 U=zeros([nx ny nz]);
 U(:,:,1:end-1)=ncread(url_uvel,'u',[i1 j1 1 it1],[nx ny nz-1 1]);
 U(:,:,end)=U(:,:,end-1);
 U=flipdim(U,3);
 
 u_roms=gridinterp_3d(XX,YY,z,U,in0,iclosest,xbu2,ybu2,zu);

 % v:
 disp('v...');
 V=zeros([nx ny nz]);
 V(:,:,1:end-1)=ncread(url_vvel,'v',[i1 j1 1 it1],[nx ny nz-1 1]);
 V(:,:,end)=V(:,:,end-1);
 V=flipdim(V,3);
 v_roms=gridinterp_3d(XX,YY,z,V,in0,iclosest,xbv2,ybv2,zv);

 % T:
 disp('T...');
 Thycom=zeros([nx ny nz]);
 Thycom(:,:,1:end-1)=ncread(url_temp,'temperature',[i1 j1 1 it1],[nx ny nz-1 1]);
 Thycom(:,:,end)=Thycom(:,:,end-1);
 Thycom=flipdim(Thycom,3);
 T_roms=gridinterp_3d(XX,YY,z,Thycom,in0,iclosest,xbr2,ybr2,zr);

 % S:
 disp('S...');
 Shycom=zeros([nx ny nz]);
 Shycom(:,:,1:end-1)=ncread(url_salt,'salinity',[i1 j1 1 it1],[nx ny nz-1 1]);
 Shycom(:,:,end)=Shycom(:,:,end-1);
 Shycom=flipdim(Shycom,3);
 S_roms=gridinterp_3d(XX,YY,z,Shycom,in0,iclosest,xbr2,ybr2,zr);
