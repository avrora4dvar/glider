%% Create obslist files
clc; clear all; 

%time 
tOut=[datenum('2011-04-01'):3:datenum('2011-10-01')]; %times for which obslist files  
tOut=[2380:3:2422]+datenum('2005-01-01'); 
tWin=3; %length window in days
dateref=datenum('2005-01-01'); 
grdfile='/home/aruba/vol2/ipasmans/Exp/Prm/grd_ow2km_r13lin_mix.nc'; 
inDir='/home/aruba/vol2/ipasmans/Exp/Obs/Exp37'; 
outDir='/home/aruba/vol2/ipasmans/Exp/Obs/Exp42'; 

%% read and select

if ~exist(outDir,'dir'); mkdir(outDir); end

con=dir(inDir); 
for k=1:length(con)
   day=regexp(con(k).name,'obslist(\d+).nc','tokens'); 
   if isempty(day); continue; end
   day=str2double(day{1}) ; 
   
   obs=read_obslist(fullfile(inDir,con(k).name)); 
   in=obs.type==6 | obs.type==7 | obs.type==8;
   if ~any(in); continue; end; 
   
   obs.lon=obs.lon(in); obs.lat=obs.lat(in); obs.z=obs.z(in); obs.t=obs.t(in); 
   obs.dir=obs.dir(in); 
   obs.type=obs.type(in); obs.sig=obs.sig(in); obs.val=obs.val(in); 
   
   write_obslist(fullfile(outDir,con(k).name),obs); 
   
   
end
