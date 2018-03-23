%% Calculate array modes
clear all; clc; 
addpath('/home/server/student/homes/ipasmans/Matlab/projects/obslist'); 

binDir='/home/aruba/vol2/ipasmans/Bin_Exp35'; 
romsDir='/home/aruba/vol2/ipasmans/Exp/Exp35/Exp37_2395'; 
outFile='/home/server/student/homes/ipasmans/Data/article_glider/ritz_vectors_exp37.mat'; 
figFile='/home/server/student/homes/ipasmans/Figures/article_glider/ritz/exp37_surface'; 

%% Grid

grdFile='/home/aruba/vol2/ipasmans/Exp/Prm/grd_ow2km_r13lin_mix.nc';
grd.lon=ncread(grdFile,'lon_rho'); 
grd.lat=ncread(grdFile,'lat_rho'); 
grd.lonu=ncread(grdFile,'lon_u'); 
grd.lonv=ncread(grdFile,'lon_v');
grd.latu=ncread(grdFile,'lat_u'); 
grd.latv=ncread(grdFile,'lat_v');
grd.zr=ncread(grdFile,'z0_r'); 
grd.h=ncread(grdFile,'h'); 

%% Find ritz vectors

r=ncread(fullfile(romsDir,'Iter1','pcg.nc'),'r');
Br=ncread(fullfile(romsDir,'Iter1','pcg.nc'),'Br'); 

[pcg,iter,obs]=iter_data(fullfile(romsDir,'Iter1','pcg.nc')); 

%% oblist

obsFile='/home/aruba/vol2/ipasmans/Exp/Obs/Exp36/obslist2395.nc'; 
obs1=read_obslist(obsFile); 

%% Construct the Ritz vectors


ritzV=iter{end}.lanczos_vector;
ritzS=iter{end}.ritz;

r=r(:,1:length(ritzS)); 
Br=Br(:,1:length(ritzS)); 

for k=1:size(ritzV,2)
    eigV(:,k)=r*ritzV(:,k); 
    eigBV(:,k)=Br*ritzV(:,k); 
    eigNorm(k)=sqrt(eigV(:,k)'*eigBV(:,k)); 
    eigBV(:,k)=eigBV(:,k).*obs.sig/eigNorm(k); 
    eigV(:,k)=eigV(:,k)/eigNorm(k); 
end

%Sort by magnitude
[~,isort]=sort(ritzS,'descend'); 

%Remove observational part
display(min(ritzS))
display(max(ritzS))
ritzS=ritzS-1; 


%% Convert from dual to primal

nVec=56; 
%isort=[1:nVec]; 
outDir=fullfile(romsDir,'Ritz'); 

vec=struct('temp',[],'salt',[],'v',[],'u',[],'zeta',[]); 
for iv=1:nVec
    display(iv); 
    iv1=isort(iv); 
    
    %Surface structure
    surface.temp(:,:,iv)=ritzS(iv1)*squeeze(...
        ncread(fullfile(outDir,sprintf('Member_%-.3d/Cv.nc',iv1)),'temp',[1,1,40,1],[Inf,Inf,1,1]) ); 
    surface.salt(:,:,iv)=ritzS(iv1)*squeeze(...
        ncread(fullfile(outDir,sprintf('Member_%-.3d/Cv.nc',iv1)),'salt',[1,1,40,1],[Inf,Inf,1,1]) ); 
    surface.u(:,:,iv)=ritzS(iv1)*squeeze(...
        ncread(fullfile(outDir,sprintf('Member_%-.3d/Cv.nc',iv1)),'u',[1,1,40,1],[Inf,Inf,1,1]) ); 
    surface.v(:,:,iv)=ritzS(iv1)*squeeze(...
        ncread(fullfile(outDir,sprintf('Member_%-.3d/Cv.nc',iv1)),'v',[1,1,40,1],[Inf,Inf,1,1]) ); 
    surface.zeta(:,:,iv)=ritzS(iv1)*squeeze(...
        ncread(fullfile(outDir,sprintf('Member_%-.3d/Cv.nc',iv1)),'zeta',[1,1,1],[Inf,Inf,1]) ); 
    
    %Find maximum in surface Ritz
    [ix1,iy1]=find(abs(surface.temp(:,:,iv))==max(reshape(abs(surface.temp(:,:,iv)),[],1)),1);
    cross.ilat(iv)=iy1; 
    
    %Cross section along maximum
    cross.lat(iv)=grd.lat(1,cross.ilat(iv)); 
    cross.zr=squeeze(grd.zr(:,cross.ilat(iv),:)); 
    cross.temp(:,:,iv)=ritzS(iv1)*squeeze(...
        ncread(fullfile(outDir,sprintf('Member_%-.3d/Cv.nc',iv1)),'temp',[1,cross.ilat(iv),1,1],[Inf,1,Inf,1]) ); 
    cross.salt(:,:,iv)=ritzS(iv1)*squeeze(...
        ncread(fullfile(outDir,sprintf('Member_%-.3d/Cv.nc',iv1)),'salt',[1,cross.ilat(iv),1,1],[Inf,1,Inf,1]) ); 
    cross.u(:,:,iv)=ritzS(iv1)*squeeze(...
        ncread(fullfile(outDir,sprintf('Member_%-.3d/Cv.nc',iv1)),'u',[1,cross.ilat(iv),1,1],[Inf,1,Inf,1]) ); 
    cross.v(:,:,iv)=.5*ritzS(iv1)*squeeze(...
        ncread(fullfile(outDir,sprintf('Member_%-.3d/Cv.nc',iv1)),'v',[1,cross.ilat(iv)-1,1,1],[Inf,1,Inf,1]) )+...
        .5*ritzS(iv1)*squeeze(...
        ncread(fullfile(outDir,sprintf('Member_%-.3d/Cv.nc',iv1)),'v',[1,cross.ilat(iv),1,1],[Inf,1,Inf,1]) ); 
    
    [min(reshape(surface.temp(:,:,iv),[],1)),max(reshape(surface.temp(:,:,iv),[],1))]
    [min(reshape(cross.temp(:,:,iv),[],1)),max(reshape(cross.temp(:,:,iv),[],1))]
        
end

%% Determine variance

for iv=1:nVec
   d=(obs.val-obs.for(:,1))./obs.sig; %Normalized innvoation vector
   ritzIn(iv)=(eigV(:,iv)'*d)/sqrt(eigV(:,iv)'*eigV(:,iv)); 
end

ritzIn=ritzIn(isort).^2;
ritzS=ritzS(isort);
p_ritzS=ritzS/sum(ritzS); 
p_ritzIn=ritzIn/sum(ritzIn); 


%% Save

save(outFile,'surface','cross','ritzS','p_ritzS','p_ritzIn','grd'); 

%% Plot 

if false
    figure();
    load(outFile);
    for iv=1:length(ritzS)
        val1=surface.temp(:,:,iv);
        val1=val1/max(abs(val1(:)));
        
        clf; hold on;
        pcolor(grd.lon,grd.lat,val1); shading flat;
        contour(grd.lon,grd.lat,grd.h,[3,200,1e3,2e3],'color',[.7,.7,.7]);
        
        in=obs1.type==6;
        plot(obs1.lon(in),obs1.lat(in),'k.');
        title(sprintf('%.1f - %.1f',p_ritzS(iv)*100,p_ritzIn(iv)*100));
        
        set(gca,'clim',[-1,1],'xlim',[-130,-123.5],'ylim',[41,50]);
        colorbar();
        
        print(sprintf('%s_%-.2d',figFile,iv),'-dpng','-r250');
        
    end
    close all;
    
end



