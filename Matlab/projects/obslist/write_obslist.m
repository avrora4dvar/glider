function [] = write_obslist(obsFile,obs,varargin)
%WRITE_OBSLIST

if exist(obsFile,'file'); delete(obsFile);end 
if ~exist(obsFile,'file')
    nccreate(obsFile,'lon','dimensions',{'K',Inf},'format','classic'); 
    nccreate(obsFile,'lat','dimensions',{'K',Inf},'format','classic'); 
    nccreate(obsFile,'z','dimensions',{'K',Inf},'format','classic'); 
    nccreate(obsFile,'time','dimensions',{'K',Inf},'format','classic'); 
    nccreate(obsFile,'type','dimensions',{'K',Inf},'format','classic','datatype','int32'); 
    nccreate(obsFile,'obs','dimensions',{'K',Inf},'format','classic'); 
    nccreate(obsFile,'sig_d','dimensions',{'K',Inf},'format','classic');
    nccreate(obsFile,'dir','dimensions',{'K',Inf},'format','classic'); 
end

if isempty(varargin)
    in=ones(size(obs.val))==1; 
else
    in=varargin{1}; 
end

ncwrite(obsFile,'lon',obs.lon(in)); 
ncwrite(obsFile,'lat',obs.lat(in)); 
ncwrite(obsFile,'z',obs.z(in)); 
ncwrite(obsFile,'time',obs.t(in)); 
ncwrite(obsFile,'type',obs.type(in)); 
ncwrite(obsFile,'obs',obs.val(in)); 
ncwrite(obsFile,'sig_d',obs.sig(in)); 
if isfield(obs,'dir')
ncwrite(obsFile,'dir',obs.dir(in)); 
end


end

