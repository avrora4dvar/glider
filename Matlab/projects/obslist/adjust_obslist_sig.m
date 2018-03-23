clear all; clc; 
addpath('/home/server/student/homes/ipasmans/Matlab/projects/obslist/'); 
obsDir='/home/aruba/vol2/ipasmans/Exp/Obs/Exp37/'; 
outDir='/home/aruba/vol2/ipasmans/Exp/Obs/Exp37/new'; 


if ~exist(outDir,'dir'); mkdir(outDir); end

%read obslist
winList=[2380:3:2422];
obsAll=struct('type',[],'sig',[],'z',[]); 
for iWin=1:length(winList)
    obs{iWin}=read_obslist(fullfile(obsDir,sprintf('obslist%-.4d.nc',winList(iWin))));
    obsAll.type=[obsAll.type;obs{iWin}.type]; 
    obsAll.sig=[obsAll.sig;obs{iWin}.sig]; 
    obsAll.z=[obsAll.z;obs{iWin}.z]; 
end

%Fields
fieldRange=[4,6,7,8,1000;...
    4,6,7,8,9999];

%Find medians
for iF=1:size(fieldRange,2)
   in=obsAll.type>=fieldRange(1,iF) & obsAll.type<=fieldRange(2,iF) & obsAll.z<6; 
   if ~any(in)
       midF(iF)=NaN;
       stdF(iF)=NaN; 
   else
       midF(iF)=median(obsAll.sig(in)); 
       stdF(iF)=std(obsAll.sig(in)); 
   end
end

%Scale. Compare var_tot-var_obs_prior with values covariance with sig=.9C
%and reference point (203,132)
model_sig=[.9,.9,.144,.0751,.0175];
stat_sig=[.9226,1.03,.5088,.1456,.0508]; %exp35_2380
%stat_sig=[.7623,.8214,.7774,.1225,.035]; %exp35_2383

%New median
for iF=1:size(fieldRange,2)
    scale(iF)=stat_sig(iF)/model_sig(iF); 
end
scale

%Remove outliers
for iObs=1:length(obs)
    for iF=1:size(fieldRange,2)
        if fieldRange(1,iF)>=6 && fieldRange(1,iF)<=7; continue; end
        in=obs{iObs}.type>=fieldRange(1,iF) & obs{iObs}.type<=fieldRange(2,iF);
        tmp=sum(in); 
        in=in & obs{iObs}.sig>midF(iF)+2*stdF(iF); 
        obs{iObs}.sig(in)=NaN; 
    end
end

for iObs=1:length(obs)
    for iF=1:size(fieldRange,2)
        if fieldRange(1,iF)>=6 && fieldRange(1,iF)<=7
            in=obs{iObs}.type==fieldRange(1,iF); 
            obs{iObs}.sig(in)=obs{iObs}.sig(in)/scale(iF); 
            obs{iObs}.sig(in)=max(...
                exp(-obs{iObs}.z(in)/100)*midF(iF)/scale(iF),...
                obs{iObs}.sig(in) ); 
        else
            in=obs{iObs}.type>=fieldRange(1,iF) & obs{iObs}.type<=fieldRange(2,iF); 
            obs{iObs}.sig(in & ~isnan(obs{iObs}.sig))=midF(iF)/scale(iF);            
        end
    end
    
    %All
    in=~isnan(obs{iObs}.sig) & ~isinf(obs{iObs}.sig); 
    obs1.lon=obs{iObs}.lon(in); 
    obs1.lat=obs{iObs}.lat(in);
    obs1.z=obs{iObs}.z(in);
    obs1.t=obs{iObs}.t(in);
    obs1.val=obs{iObs}.val(in);
    obs1.sig=obs{iObs}.sig(in);
    obs1.dir=obs{iObs}.dir(in);
    obs1.type=obs{iObs}.type(in);
    write_obslist(fullfile(outDir,sprintf('obslist%-.4d.nc',winList(iObs))),obs1); 
    
  
end






