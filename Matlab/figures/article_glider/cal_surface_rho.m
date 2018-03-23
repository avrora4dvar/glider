%% Read surfaces
clear all; clc; 
addpath('/home/server/student/homes/ipasmans/Matlab/roms/seawater'); 

%roms directories
romsDir={'/home/aruba/vol2/ipasmans/Exp/Exp31/Exp35_ana/',...
   '/home/aruba/vol2/ipasmans/Exp/Exp31/Exp36_ana/',...
   '/home/aruba/vol2/ipasmans/Exp/Exp31/Exp37_ana/',...
   '/home/aruba/vol2/ipasmans/Exp/Exp31/Exp40_ana/',...
   };

tWin=[2392:3:2410]; 
dateRef=datenum('2005-01-01'); 

outFile='/home/server/student/homes/ipasmans/Data/article_glider/surface_rho.mat'; 

%% Read temp,salt,rho

model=[];
for iExp=1:length(romsDir)
    
    model1=struct('dir',romsDir{iExp},'t',[],'temp',[],'salt',[],'rho',[]); 
    for iWin=1:length(tWin)
        fname=fullfile(romsDir{iExp},sprintf('ocean_avg_%-.4d.nc',tWin(iWin)));
        
        t1=ncread(fname,'ocean_time')/24/3600+dateRef; 
        model1.t=[model1.t,t1(:)']; 
        temp1=ncread(fname,'temp',[1,1,40,1],[Inf,Inf,1,Inf]); 
        model1.temp=cat(3,model1.temp,squeeze(temp1)); 
        salt1=ncread(fname,'salt',[1,1,40,1],[Inf,Inf,1,Inf]); 
        model1.salt=cat(3,model1.salt,squeeze(salt1)); 
        val1=sw_pden(salt1(:),temp1(:),0*salt1(:),0); 
        model1.rho=cat(3,model1.rho,reshape(val1,[size(salt1,1),size(salt1,2),size(salt1,4)])); 
                
    end
    
    model=[model,model1]; 
end

save(outFile,'model'); 
