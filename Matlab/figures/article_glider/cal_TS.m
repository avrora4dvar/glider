%% Read temperature and salinity observations at NH10 and compare with models
clear all; clc; 

addpath('..'); 
addpath('/home/server/student/homes/ipasmans/Matlab/roms/seawater/'); 
addpath('/home/server/student/homes/ipasmans/Matlab/figures/io/'); 

%directory with glider files
obsDir='/home/server/student/homes/ipasmans/Obs/glider/osu/';

%grid file
grid_file='/home/aruba/vol2/ipasmans/Exp/Prm/grd_ow2km_r13lin_mix.nc';

romsDir={'/home/aruba/vol2/ipasmans/Exp/Exp35/Exp36_2392/NL',...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp36a_2392/NL',...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp40_ana'}; 

romsDir={'/home/aruba/vol2/ipasmans/Exp/Exp35/Exp35_ana',...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp36_ana',...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp37_ana',...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp38_ana',...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp39_ana',...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp40_ana'}; 


    

%time boundaries
tLim=[2392,2413]+datenum('2005-01-01');
dLim=[0,250]; 

%output
outName='/home/server/student/homes/ipasmans/Data/article_glider/TS_glider_ana.mat'; 

%% read grid

grd.lon=ncread(grid_file,'lon_rho'); 
grd.lat=ncread(grid_file,'lat_rho'); 
grd.mask=ncread(grid_file,'mask_rho'); 


%% read observations

obs=io_glider_osu(obsDir,'dir'); 
obs.source='http://gliderfs2.coas.oregonstate.edu/gliderweb/archive/gridded/slocum_gridded_data.php'; 
in=obs.t>=tLim(1) & obs.t<=tLim(2); 
obs.lon=obs.lon(in); obs.lat=obs.lat(in); obs.t=obs.t(in); 
obs.temp=obs.temp(in,:); obs.salt=obs.salt(in,:); obs.p=obs.p(in,:); obs.depth=obs.depth(in,:); 
obs.sig_temp=obs.sig_temp(in,:); obs.sig_salt=obs.sig_salt(in,:); obs.sig_p=obs.sig_p(in,:); 

loc=struct('lon',repmat(obs.lon,[1,length(obs.cell_depth)]),...
    'lat',repmat(obs.lat,[1,length(obs.cell_depth)]),...
    'depth',obs.depth,'t',repmat(obs.t,[1,length(obs.cell_depth)])); 
loc.lon=reshape(loc.lon,1,[]); 
loc.lat=reshape(loc.lat,1,[]); 
loc.depth=reshape(loc.depth,1,[]); 
loc.t=reshape(loc.t,1,[]); 

if length(loc.t)==0
    error('No observations found'); 
end
%% interpolate models

model=[]; 
for iMod=1:length(romsDir)
    
    
    display(sprintf('Interpolating temp from %s',romsDir{iMod})); 
    mod1=struct('temp',[],'salt',[]); 
    out1=get_roms_3d(grid_file,romsDir{iMod},'temp',grd,loc);
    temp1=reshape(out1.val,size(obs.temp)); 
    display(sprintf('Interpolating salt from %s',romsDir{iMod})); 
    out1=get_roms_3d(grid_file,romsDir{iMod},'salt',grd,loc);
    salt1=reshape(out1.val,size(obs.salt)); 
    z1=reshape(out1.z,size(obs.depth)) ; 
    
    model=[model,struct(...
        'temp',temp1,'salt',salt1,...
         'romsDir',romsDir{iMod},'z',z1 ) ];
end




%% save

display('Saving output'); 
out_file=sprintf('%s',outName); 
save(out_file,'loc','obs','model'); 