%% Calculate surface currents and vorticity surface layer
clear all; clc;
addpath('../'); 

grdFile='/home/aruba/vol2/ipasmans/Exp/Prm/grd_ow2km_r13lin_mix.nc'; 
romsDir={...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp36_ana/',...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp37_ana/'}; 
winList=[2392:3:2395]; 
dateRef=datenum('2005-01-01'); 
outFile='/home/server/student/homes/ipasmans/Data/article_glider/exp35_vort_ana.mat'; 

%% Read grid

grd.lon=ncread(grdFile,'lon_rho'); 
grd.lat=ncread(grdFile,'lat_rho'); 
grd.h=ncread(grdFile,'h'); 
grd.f=sin(deg2rad(grd.lat))*2*7.2921e-5; 

%% Read fields

model=[];
for iExp=1:length(romsDir)
    model1.dir=romsDir{iExp};
    model1.u=[]; model1.v=[]; model1.rvort=[]; model1.t=[]; 
  
    for iWin=1:length(winList)
        fname=fullfile(romsDir{iExp},sprintf('ocean_avg_%-.4d.nc',winList(iWin))); 
        display(fname); 
        u=ncread(fname,'u',[1,1,40,1],[Inf,Inf,1,Inf]); 
        v=ncread(fname,'v',[1,1,40,1],[Inf,Inf,1,Inf]); 
        t1=ncread(fname,'ocean_time')/3600/24+dateRef;
       
        for it=1:length(t1)
            model1.t=[model1.t,t1(it)]; 
            
            val=nan(size(grd.lon,1),size(grd.lon,2)); 
            val=squeeze(.5*u(1:end-1,2:end-1,1,it)+.5*u(2:end,2:end-1,1,it)); 
            model1.u=cat(3,model1.u,val); 
            val=squeeze(.5*v(2:end-1,1:end-1,1,it)+.5*v(2:end-1,2:end,1,it)); 
            model1.v=cat(3,model1.v,val); 
            
            val=uv2vort(grdFile,squeeze(u(:,:,1,it)),squeeze(v(:,:,1,it))); 
            model1.rvort=cat(3,model1.rvort,val); 
        end
        
    end
    model=[model,model1]; 
end

%% Save

save(outFile,'model','grd'); 



  