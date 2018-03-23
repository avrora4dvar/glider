%% Calculate rho surface
addpath('/home/server/student/homes/ipasmans/Matlab/roms/seawater/'); 
addpath('/home/server/student/homes/ipasmans/Matlab/figures'); 


%Input
dateRef=datenum('2005-01-01'); 
expDir={'/home/aruba/vol2/ipasmans/Exp/Exp35/Exp35_for',...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp36_for/',...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp37_for',...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp38_for',...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp40_ana'}; 
grdFile='/home/aruba/vol2/ipasmans/Exp/Prm/grd_ow2km_r13lin_mix.nc'; 

%Output
t=[2392.5:2412.5]+dateRef; 
contours=[1023.5:1027.5];
outFile='/home/server/student/homes/ipasmans/Data/article_glider/rho_contours_for.mat'; 


%% Read grid

grd.sigma=[-39.5:-.5]/40; 
grd.sigmaw=[-40:0]/40; 
grd.lon=ncread(grdFile,'lon_rho'); 
grd.lat=ncread(grdFile,'lat_rho'); 
grd.mask=ncread(grdFile,'mask_rho'); 
grd.h=ncread(grdFile,'h'); 
grd.s=size(grd.lon); 

% grd3.sigma=repmat(reshape(grd.sigma,1,1,[]),[size(grd.lon,1),size(grd.lon,2),1]); 
% grd3.lon=repmat(grd.lon,[1,1,size(grd3.sigma,3)]); 
% grd3.lat=repmat(grd.lat,[1,1,size(grd3.sigma,3)]); 
% grd3.lon=permute(grd3.lon,[2,1,3]); 
% grd3.lat=permute(grd3.lat,[2,1,3]); 
% grd3.sigma=permute(grd3.sigma,[2,1,3]); 

%% Read date

%Output
model=[]; 

matlabpool('open',6); 
for iMod=1:length(expDir)
    con=dir(expDir{iMod}); 
    
    %Output
    model1.lon=grd.lon; 
    model1.lat=grd.lat; 
    model1.contours=contours; 
    model1.t=t; 
    model1.dir=expDir{iMod}; 
    model1.z=nan(size(grd.lon,1),size(grd.lon,2),length(contours),length(t)); 
     
    for iFile=1:length(con)
        %Read times from file
        if isempty(regexp(con(iFile).name,'ocean_avg_(\d+).nc','once')); continue; end
        fname=fullfile(expDir{iMod},con(iFile).name); 
        display(fname); 
        t1=ncread(fname,'ocean_time')/24/3600+dateRef; 
        
        for itIn=1:length(t1)
            %Find time in output
           itOut=find(t==t1(itIn),1); 
           if isempty(itOut); continue; end
           display(datestr(t1(itIn),'yyyy-mm-dd HH:MM')); 
           
           %Read fields
           temp1=squeeze(ncread(fname,'temp',[1,1,1,itIn],[Inf,Inf,Inf,1])); 
           salt1=squeeze(ncread(fname,'salt',[1,1,1,itIn],[Inf,Inf,Inf,1])); 
           zeta1=squeeze(ncread(fname,'zeta',[1,1,itIn],[Inf,Inf,1])); 
           z1=squeeze(roms_sigma2z(grd.sigma,grd.h,zeta1)); 
           z2=squeeze(roms_sigma2z(grd.sigmaw,grd.h,zeta1)); 
           
           %Calculate density
           rho1=nan(size(temp1)); 
           for i3=1:size(temp1,3)
               p1=sw_pres(-squeeze(z1(:,:,i3)),grd.lat); 
               rho1(:,:,i3)=sw_pden(...
                   squeeze(salt1(:,:,i3)),squeeze(temp1(:,:,i3)),...
                   squeeze(p1),0); 
           end
           
           %Stable
           rho1=stable_rho(rho1,diff(z2,1,3),'min_step',1e-6); 
           
           %Find z-coordinate contours
           ztmp=nan(length(contours),numel(grd.lon)); 
           temptmp=nan(length(contours),numel(grd.lon)); 
           salttmp=nan(length(contours),numel(grd.lon)); 
           parfor i0=1:numel(grd.lon)
               [i1,i2]=ind2sub(grd.s,i0); 
               if grd.mask(i0)==0; continue; end
               
               ztmp(:,i0)=interp1(squeeze(rho1(i1,i2,:)),squeeze(z1(i1,i2,:)),contours);
                end
           model1.z(:,:,:,itOut)=permute(reshape(ztmp,[],grd.s(1),grd.s(2)),[2,3,1]); 
         
        end %itIn
    end %iFile
    
    %Reshape output
    model=[model,model1]; 

end %iMod
matlabpool('close'); 

%Save
save(outFile,'model'); 
