function [obs] = read_obslist(obsFile,varargin)
%READ_OBSLIST [obslist_struct]=read_obslist(obslist_file,<window_directory>)

info=ncinfo(obsFile); 


%% read obslist

obs.lon=ncread(obsFile,'lon'); 
obs.lat=ncread(obsFile,'lat');
obs.z=ncread(obsFile,'z'); 
obs.t=ncread(obsFile,'time'); 
obs.type=ncread(obsFile,'type'); 
obs.val=ncread(obsFile,'obs'); 
obs.sig=ncread(obsFile,'sig_d'); 
if ismember('dir',[info.Variables(:).Name])
    obs.dir=ncread(obsFile,'dir'); 
end

ref=regexp(obsFile,'obslist(\d+).nc','tokens'); 
if ~isempty(ref)
    obs.refDate=str2double(ref{1}{1}); 
    obs.datenum=datenum('2005-01-01')+obs.refDate+obs.t/24/3600;
end

%% Read samples


for iFile=1:2:length(varargin)
    romsDir=fullfile(varargin{iFile},sprintf('%s_%-.4d',varargin{iFile+1},obs.refDate));
    obs.for=ncread(fullfile(romsDir,'Iter0','sample_for.nc'),'K');
    obs.ana=ncread(fullfile(romsDir,'Iter1','sample_ana.nc'),'K');
    
    %Try reading ensemble
    for k=1:999
        memberDir=fullfile(romsDir,'Ens',sprintf('Member_%-.3d',k)); 
        if exist(fullfile(memberDir,'sample_for.nc'),'file')
            obs.ensFor(:,k)=ncread(fullfile(memberDir,'sample_for.nc'),'K');
        end
        if exist(fullfile(memberDir,'sample_ana.nc'),'file')
            obs.ensAna(:,k)=ncread(fullfile(memberDir,'sample_ana.nc'),'K');
        end 
    end
end




end

