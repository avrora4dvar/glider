%% create discharge time series for Columbia River

%% Input
tOut=[datenum('2010-12-31'):datenum('2012-01-01')]+.5;
file_out='ow2km_rivers_2011.nc'; 
nRivers=2; nLayers=40; 
dateIni=datenum('2005-01-01'); 

%% Create output file 

if exist(file_out,'file'); delete(file_out); end

nccreate(file_out,'river','Format','classic','Datatype','int32','Dimensions',{'river',nRivers}); 
ncwriteatt(file_out,'river','long_name','river runoff identification number'); 
ncwriteatt(file_out,'river','units','nondimensional');
ncwriteatt(file_out,'river','field','river,scalar'); 

nccreate(file_out,'river_Xposition','Format','classic','Datatype','int32','Dimensions',{'river',nRivers}); 
ncwriteatt(file_out,'river_Xposition','long_name','river XI-position at RHO-points'); 
ncwriteatt(file_out,'river_Xposition','units','nondimensional');
ncwriteatt(file_out,'river_Xposition','field','river,scalar'); 

nccreate(file_out,'river_Eposition','Format','classic','Datatype','int32','Dimensions',{'river',nRivers}); 
ncwriteatt(file_out,'river_Eposition','long_name','river ETA-position at RHO-points'); 
ncwriteatt(file_out,'river_Eposition','units','nondimensional');
ncwriteatt(file_out,'river_Eposition','field','river,scalar'); 

nccreate(file_out,'river_direction','Format','classic','Datatype','int32','Dimensions',{'river',nRivers}); 
ncwriteatt(file_out,'river_direction','long_name','river runoff direction'); 
ncwriteatt(file_out,'river_direction','units','nondimensional');
ncwriteatt(file_out,'river_direction','field','river,scalar'); 

%0=all tracers are off,1=only temp,only salinity,3=temp and salinity,
nccreate(file_out,'river_flag','Format','classic','Datatype','int32','Dimensions',{'river',nRivers}); 
ncwriteatt(file_out,'river_flag','long_name','river runoff tracer flag'); 
ncwriteatt(file_out,'river_flag','units','nondimensional');
ncwriteatt(file_out,'river_flag','field','river,scalar'); 

nccreate(file_out,'river_Vshape','Format','classic','Datatype','single','Dimensions',{'river',nRivers,'s_rho',nLayers}); 
ncwriteatt(file_out,'river_Vshape','long_name','river runoff mass transport vertical profile'); 
ncwriteatt(file_out,'river_Vshape','units','nondimensional');
ncwriteatt(file_out,'river_Vshape','field','river,scalar'); 

nccreate(file_out,'river_time','Format','classic','Datatype','double','Dimensions',{'time',Inf}); 
ncwriteatt(file_out,'river_time','long_name','river runoff time'); 
ncwriteatt(file_out,'river_time','units','days since 2005-01-01 00:00:00');

nccreate(file_out,'river_transport','Format','classic','Datatype','single','Dimensions',{'river',nRivers,'time',Inf}); 
ncwriteatt(file_out,'river_transport','long_name','river runoff vertically integrated mass transport'); 
ncwriteatt(file_out,'river_transport','units','meter3 second-1');
ncwriteatt(file_out,'river_transport','field','river_transport, scalar, series'); 
ncwriteatt(file_out,'river_transport','time','river_time'); 

nccreate(file_out,'river_temp','Format','classic','Datatype','single','Dimensions',{'river',nRivers,'s_rho',nLayers,'time',Inf}); 
ncwriteatt(file_out,'river_temp','long_name','river runoff potential temperature'); 
ncwriteatt(file_out,'river_temp','units','Celsius');
ncwriteatt(file_out,'river_temp','field','river_temp, scalar, series'); 
ncwriteatt(file_out,'river_temp','time','river_time'); 

nccreate(file_out,'river_salt','Format','classic','Datatype','single','Dimensions',{'river',nRivers,'s_rho',nLayers,'time',Inf}); 
ncwriteatt(file_out,'river_salt','long_name','river runoff salinity'); 
ncwriteatt(file_out,'river_salt','units','ppt');
ncwriteatt(file_out,'river_salt','field','river_salt, scalar, series'); 
ncwriteatt(file_out,'river_salt','time','river_time'); 



%% Read discharge data Columbia (USGS 14246900)

%storage
dis=struct('t',[],'val',[]); 

%read file
fid=fopen('CR_2011_discharge.txt'); 
while ~feof(fid)
    line=fgetl(fid); 
    line=regexp(line,'\t','split');
    
    if isempty(line) || ~strcmp(line{1},'USGS')
        continue;
    end
    
    dis.t=[dis.t,datenum(line{3},'yyyy-mm-dd')+.5];
    dis.val=[dis.val,str2num(line{4})]; 
end
fclose(fid); 
Columbia River
%convert time to GMT
dis.t=dis.t+7/24; 
%convert m3 s-1
dis.val=dis.val/35.3147248277; 

calculate daily values
dis.out=interp1(dis.t,dis.val,tOut); 

%% Read temperature Columbia (CMOP Woody Island)

%storage
temp=struct('t',[],'val',[]); 

%read file
fid=fopen('CR_2011_temp.csv'); 
line=fgetl(fid); 
while ~feof(fid)
    line=fgetl(fid); 
    line=regexp(line,',','split');
    
    if isempty(line) 
        continue;
    end
    
     temp.t=[temp.t,datenum(line{1})];
     temp.val=[temp.val,str2num(line{2})]; 
end
fclose(fid); 

%order ascending in time
tmp=[temp.t(:),temp.val(:)]; 
tmp=sortrows(tmp); 
temp.t=tmp(:,1)'; temp.val=tmp(:,2)'; 

%convert to GMT
temp.t=temp.t+7/24; 

%calculate daily values
for it=1:length(tOut)
    inPeriod=temp.t>=tOut(it)-.5&temp.t<tOut(it)+.5;
    temp.out(it)=trapz(temp.t(inPeriod),temp.val(inPeriod))/trapz(temp.t(inPeriod),ones(1,sum(inPeriod))); 
end

%% Columbia River

%location
ncwrite(file_out,'river',[1,2],[1]); 
ncwrite(file_out,'river_Xposition',[272,272],[1]); 
ncwrite(file_out,'river_Eposition',[311,312],[1]); 
ncwrite(file_out,'river_direction',[0,0],[1]); 
ncwrite(file_out,'river_flag',[3,3],[1]); 
%vertical distribution
val=ones(1,nLayers)/nLayers;
ncwrite(file_out,'river_Vshape',val,[1,1]);
ncwrite(file_out,'river_Vshape',val,[2,1]);
%time
ncwrite(file_out,'river_time',tOut-dateIni,[1]); 
%transport
ncwrite(file_out,'river_transport',-.5*dis.out,[1,1]); 
ncwrite(file_out,'river_transport',-.5*dis.out,[2,1]); 
%temperature
val=repmat(temp.out,[1,nLayers,1]); 
ncwrite(file_out,'river_temp',val,[1,1,1]); 
ncwrite(file_out,'river_temp',val,[2,1,1]); 
%salinity
val=.3*ones(1,nLayers,length(tOut)); 
ncwrite(file_out,'river_salt',val,[1,1,1]); 
ncwrite(file_out,'river_salt',val,[2,1,1]); 