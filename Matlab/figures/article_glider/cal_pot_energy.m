%% Density
clear all; clc; 
addpath('..'); 

romsDir={...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp35_ana/',...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp37_ana/',...
    };

%% Grid

grdFile='/home/aruba/vol2/ipasmans/Exp/Prm/grd_ow2km_r13lin_mix.nc';
grd.mask=ncread(grdFile,'mask_rho'); 
grd.lon=ncread(grdFile,'lon_rho'); 
grd.lat=ncread(grdFile,'lat_rho'); 
grd.zr=ncread(grdFile,'z0_r'); 
grd.zw=ncread(grdFile,'z0_w'); 
grd.dz=diff(grd.zw,[],3); 
grd.dz2=.5*grd.zw(:,:,2:end).^2-.5*grd.zw(:,:,1:end-1).^2; 
grd.dz3=-1/3*grd.zw(:,:,2:end).^3+1/3*grd.zw(:,:,1:end-1).^3; 

% in=repmat(grd.mask,[1,1,41]); in(in==0)=NaN; in=grd.zw.*in;
% for i1=1:size(grd.zw,3)
%     z(i1,1)=nanmean(reshape(in(:,:,i1),[],1));
% end
% z=[min(grd.zr(:));z];
% z=interp1([1:length(z)],z,[1:.8:length(z)]);
z=[-500:4:0]; 
grd.mask(grd.lon>-123.92)=0; 
grd.mask(grd.lon>-124.669 & grd.lat>48.037)=0; 
grd.mask(grd.lon>-124.957 & grd.lat>48.759)=0; 
grd.mask3=repmat(grd.mask,[1,1,size(grd.zr,3)]); 
grd.mask3(grd.mask3==0)=NaN; 

%% 

surface=cell(length(romsDir),1); 
crossZ=cell(length(romsDir),1); 
crossM=cell(length(romsDir),1); 
t=[]; 
for iWin=[2410:3:2410]
    for iMod=1:length(romsDir)
        display([iMod,iWin]);
        fname=fullfile(romsDir{iMod},sprintf('ocean_avg_%-.4d.nc',iWin));
        t1=ncread(fname,'ocean_time')/24/3600+datenum('2005-01-01'); 
        temp1=ncread(fname,'temp'); 
        salt1=ncread(fname,'salt'); 
        s=size(salt1); 
        
        for it=1:1 %size(temp1,4)
            rho1{iMod}(:,:,:,it)=sw_pden(salt1(:,:,:,it),temp1(:,:,:,it),-grd.zr,0*grd.zr); 
        
            %Turn into stable profile
            rho1{iMod}(:,:,:,it)=stable_rho(rho1{iMod}(:,:,:,it),grd.dz,'min_step',1e-6); 
            
            rho1{iMod}(:,:,:,it)=rho1{iMod}.*grd.mask3; 
        end
        
     
        Epot=9.81*bsxfun(@times,grd.dz2,abs(rho1{iMod}));
        for it=1:1
           surface{iMod}=cat(3,surface{iMod},squeeze(nansum(Epot,3))); 
        end
        

        Epot=Epot./diff(grd.zw,[],3); 
        for i1=1:length(z)-1
           in=grd.zw(:,:,2:end)>=z(i1); 
           in=in & grd.zw(:,:,1:end-1)<z(i1); 
           Epot1=Epot; Epot1(~in)=NaN; 
           
           crossZ{iMod}=cat(2,crossZ{iMod},...
               2e3*nansum(nansum(Epot1,3),2)); 
           crossM{iMod}=cat(2,crossM{iMod},...
               permute(2e3*nansum(nansum(Epot1,3),1),[2,1]));  
        end 
        
        %Calculate pot. energy
        if iMod==1
            t=cat(1,t,t1);
        end
        
    end
    

end

%%

z=.5*z(1:end-1)+.5*z(2:end); 

%lon=grd.lon; lat=grd.lat; 
save('/home/server/student/homes/ipasmans/Data/article_glider/Epot_2410.mat','z','surface','crossM','crossZ','grd'); 
%load('/home/server/student/homes/ipasmans/Data/article_glider/Epot.mat','z','surface','crossM','crossZ','grd'); 