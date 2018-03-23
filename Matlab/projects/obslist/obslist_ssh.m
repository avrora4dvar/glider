function [out]=obslist_ssh(grdfile,dataDir,tide,outfile,tLim,mask)
% Data type format:
% fnnn, where f is the satellite, nnn track number
% f: j1      | 1
%    j2      | 2
%    en      | 3
%    cryosat | 4
%    altika  | 5

display(sprintf('%s: Starting obslist_ssh',datestr(now(),'dd mmm HH:MM')));
addpath('/home/server/student/homes/ipasmans/Matlab/roms/t_tide'); 
addpath('/home/server/student/homes/ipasmans/Matlab/figures'); 

%options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Based on difference 24h mean Exp26f_2353
%dSST/dxi: .0029 m dSST/deta;.0022
%dSST/dxi: .0017 m dSST/deta:NaN
sig0=.003;

sig_d=0.02; % obs std dev (m) (was .02)
flag_tide=0; %(0=no tide, 1=add tide)
cycle_length=20; %minimal number of observations in cycle 
assim_type='daily_average'; 

%read grid
grd.lon=ncread(grdfile,'lon_rho');
grd.lat=ncread(grdfile,'lat_rho'); 
grd.mask=ncread(grdfile,'mask_rho')==1; 

%limits (1 cell from boundary)
xLim=[grd.lon(2,1),grd.lon(end-1,end)]; 
yLim=[grd.lat(1,2),grd.lat(end,end-1)]; 


out=struct('t',[],'lon',[],'lat',[],'cycle',[],'ssh',[],'sig',[]); 

%% read jason-1 (1000)

if true

satDir=fullfile(dataDir,'j1'); 
if exist(satDir,'dir')
    conSat=dir(satDir); 
    for i0=1:length(conSat)
       cycleDir=fullfile(satDir,conSat(i0).name); 
       conCycle=dir(cycleDir); 
       
       for i1=1:length(conCycle)
          if ~isempty(regexp(conCycle(i1).name,'.nc','once'))
              fname=fullfile(cycleDir,conCycle(i1).name); 
              [bndB,bndE]=jason_tLim(fname);
              if tLim(1)>bndE || tLim(2)<bndB
                  continue;
              end
              
              display(sprintf('reading %s',fname)); 
              time=ncread(fname,'time')/24/3600+datenum('2000-01-01');
              
              %indices with time limits
              bndB=find(time>=tLim(1),1,'first'); 
              bndE=find(time<=tLim(2),1,'last'); 
              if isempty(bndB) || isempty(bndE)
                continue; 
              end
              
              %read values
              lon=ncread(fname,'lon',[bndB],[bndE-bndB+1]); lon(lon>180)=lon(lon>180)-360; 
              lat=ncread(fname,'lat',[bndB],[bndE-bndB+1]); 
              flag=ncread(fname,'alt_quality_flag',[bndB],[bndE-bndB+1])==0 &...
                  ncread(fname,'rad_quality_flag',[bndB],[bndE-bndB+1])==0 & ...
                  ncread(fname,'geophysical_quality_flag',[bndB],[bndE-bndB+1])==0 & ...
                  ncread(fname,'ecmwf_meteo_map_avail',[bndB],[bndE-bndB+1])==0 &...
                  ncread(fname,'rain_flag',[bndB],[bndE-bndB+1])==0;     %after exp48
              ssh=ncread(fname,'mean_topography',[bndB],[bndE-bndB+1])+...
                  ncread(fname,'ssha',[bndB],[bndE-bndB+1])+...
                  flag_tide*ncread(fname,'ocean_tide_sol1',[bndB],[bndE-bndB+1]); 
              time=time(bndB:bndE); 
              cycleNo=double(ncreadatt(fname,'/','cycle_number')); 
          
              %select output
              maskVal=interp2(grd.lon',grd.lat',mask',lon,lat);
              in =...
                  lon>=xLim(1) & lon<=xLim(2) & ...
                  lat>=yLim(1) & lat<=yLim(2) & ...
                  maskVal==1;
              in=in(:);
              
              %add to output
              out.t=[out.t;time(in)];
              out.lon=[out.lon;lon(in)];
              out.lat=[out.lat;lat(in)];
              out.ssh=[out.ssh;ssh(in)];
              out.cycle=[out.cycle;cycleNo*ones(sum(in),1)+1e3];
              out.sig=[out.sig;ones(sum(in),1)*.02]; 
               
          end
       end
    end
    
end

end

%% read jason-2 (2000)

if true

satDir=fullfile(dataDir,'j2'); 
if exist(satDir,'dir')
    conSat=dir(satDir); 
    for i0=1:length(conSat)
       cycleDir=fullfile(satDir,conSat(i0).name); 
       conCycle=dir(cycleDir); 
       
       for i1=1:length(conCycle)
          if ~isempty(regexp(conCycle(i1).name,'.nc','once'))
              fname=fullfile(cycleDir,conCycle(i1).name); 
              [bndB,bndE]=jason_tLim(fname); 
              if tLim(1)>bndE || tLim(2)<bndB
                  continue;
              end
              
              display(sprintf('reading %s',fname)); 
              time=ncread(fname,'time')/24/3600+datenum('2000-01-01');
              
              %indices with time limits
              bndB=find(time>=tLim(1),1,'first'); 
              bndE=find(time<=tLim(2),1,'last'); 
              if isempty(bndB) || isempty(bndE)
                continue; 
              end
              
              %read values
              lon=ncread(fname,'lon',[bndB],[bndE-bndB+1]); lon(lon>180)=lon(lon>180)-360; 
              lat=ncread(fname,'lat',[bndB],[bndE-bndB+1]); 
              flag=ncread(fname,'alt_quality_flag',[bndB],[bndE-bndB+1])==0 &...
                  ncread(fname,'rad_quality_flag',[bndB],[bndE-bndB+1])==0 & ...
                  ncread(fname,'geophysical_quality_flag',[bndB],[bndE-bndB+1])==0 & ...
                  ncread(fname,'ecmwf_meteo_map_avail',[bndB],[bndE-bndB+1])==0 &...
                  ncread(fname,'rain_flag',[bndB],[bndE-bndB+1])==0;     %after exp48
              ssh=ncread(fname,'mean_topography',[bndB],[bndE-bndB+1])+...
                  ncread(fname,'ssha',[bndB],[bndE-bndB+1])+...
                  flag_tide*ncread(fname,'ocean_tide_sol1',[bndB],[bndE-bndB+1]);
              time=time(bndB:bndE);
              cycleNo=double(ncreadatt(fname,'/','cycle_number')); 
              
              %select output
              maskVal=interp2(grd.lon',grd.lat',mask',lon,lat);
              in =...
                  lon>=xLim(1) & lon<=xLim(2) & ...
                  lat>=yLim(1) & lat<=yLim(2) & ...
                  maskVal==1;
              in=in(:);
             
              
              %add to output
              out.t=[out.t;time(in)];
              out.lon=[out.lon;lon(in)];
              out.lat=[out.lat;lat(in)];
              out.ssh=[out.ssh;ssh(in)];
              out.cycle=[out.cycle;cycleNo*ones(sum(in),1)+2e3];
              out.sig=[out.sig;ones(sum(in),1)*.02]; 
             
          end
       end
    end

end

end
%% read Cryosat (3000)

if true

satDir=fullfile(dataDir,'c2p');
if exist(satDir,'dir')
    conSat=dir(satDir);
    for i0=1:length(conSat)
        if ~isempty(regexp(conSat(i0).name,'.nc','once'))
            fname=fullfile(satDir,conSat(i0).name);
            [bndB,bndE]=cryosat_tLim(fname); 
            if tLim(1)>bndE || tLim(2)<bndB
                continue;
            end
            
            display(sprintf('reading %s',fname));
            time=ncread(fname,'time')/24/3600+datenum('1985-01-01');
            
            %indices with time limits
            bndB=find(time>=tLim(1),1,'first');
            bndE=find(time<=tLim(2),1,'last');
            if isempty(bndB) || isempty(bndE)
                continue;
            end
            
            %read values
            lon=ncread(fname,'lon',[bndB],[bndE-bndB+1]); lon(lon>180)=lon(lon>180)-360;
            lat=ncread(fname,'lat',[bndB],[bndE-bndB+1]);
            flag=ncread(fname,'flags',[bndB],[bndE-bndB+1])==0;
            ssh=ncread(fname,'ssha',[bndB],[bndE-bndB+1])+...
                ncread(fname,'mss_cnescls11',[bndB],[bndE-bndB+1])-...
                ncread(fname,'geoid_egm2008',[bndB],[bndE-bndB+1])+...
                flag_tide*ncread(fname,'tide_ocean_fes12',[bndB],[bndE-bndB+1]);
            time=time(bndB:bndE);
            cycleNo=double(ncreadatt(fname,'/','cycle_number'));
            
            
            %select output
            ssh(~flag)=NaN; 
            maskVal=interp2(grd.lon',grd.lat',mask',lon,lat);
            in =...
                lon>=xLim(1) & lon<=xLim(2) & ...
                lat>=yLim(1) & lat<=yLim(2) & ...
                maskVal==1;
            in=in(:);
            
            %add to output
            out.t=[out.t;time(in)];
            out.lon=[out.lon;lon(in)];
            out.lat=[out.lat;lat(in)];
            out.ssh=[out.ssh;ssh(in)];
            out.cycle=[out.cycle;cycleNo*ones(sum(in),1)+3e3];
            out.sig=[out.sig;ones(sum(in),1)*.02]; %.034 according to article
            
        end
    end
end

end

%% read Envisat (4000)

if true

satDir=fullfile(dataDir,'envisat');
if exist(satDir,'dir')
    conSat=dir(satDir);
    for i0=1:length(conSat)
        if ~isempty(regexp(conSat(i0).name,'.mat','once'))
            fname=fullfile(satDir,conSat(i0).name);
            obs1=load(fname); obs1=obs1.obs; 
            
            %read time
            time=obs1.t; 
            in=time>=tLim(1) & time<=tLim(2); 
            if ~any(in); continue; end
            time=time(in); 
            
            %Read data 
            display(sprintf('reading %s',fname));
            lon=obs1.lon(in); lon(lon>180)=lon(lon>180)-360;
            lat=obs1.lat(in); 
            ssh=obs1.ssha(in); %SSHA with respect to geoid
            cycleNo=obs1.cycle(in); 
            
            %select output
            maskVal=interp2(grd.lon',grd.lat',mask',lon,lat);
            in =...
                lon>=xLim(1) & lon<=xLim(2) & ...
                lat>=yLim(1) & lat<=yLim(2) & ...
                maskVal==1;
            in=in(:);
            
            %add to output
            out.t=[out.t;time(in)];
            out.lon=[out.lon;lon(in)];
            out.lat=[out.lat;lat(in)];
            out.ssh=[out.ssh;ssh(in)];
            out.cycle=[out.cycle;cycleNo(in)+4e3];
            out.sig=[out.sig;ones(sum(in),1)*.02]; %.054 arcording to article 
            
        end
    end
end

end

%% break cycles

[cUni,c2uni,uni2c]=unique(out.cycle); 
for iCycle=1:length(cUni) 
    in=uni2c==iCycle; 
    
    time=out.t(in);
    dCycle=zeros(size(time)); 
    
    iSplit=find(time>time(1)+1/24,1,'first'); 
    while ~isempty(iSplit)
        dCycle(iSplit:end)=dCycle(iSplit:end)+200;
        iSplit=find(time>time(iSplit)+1/24,1,'first'); 
    end
    
    dCycle=mod(dCycle,1e3); 
    out.cycle(in)=out.cycle(in)+dCycle; 
        
end


%% average in time

if strcmpi(assim_type,'instant')
    out.t=out.t;
elseif strcmpi(assim_type,'daily_average')
    [cUni,c2uni,uni2c]=unique(out.cycle); 
    for iCycle=1:length(cUni)
       in=uni2c==iCycle; 
       %out.t(in)=round(mean(out.t(in))*24)/24; 
       out.t(in)=mean(out.t(in)); 
    end
else
    error('unknown assim_type'); 
end

%% Remove outliers

[cUni,c2uni,uni2c]=unique(out.cycle);
for iCycle=1:length(cUni)
    in=uni2c==iCycle; 
           
    %Remove outliers
    [out.ssh(in),out.p(in)]=outlier_bayes(out.sig(in),.5,out.ssh(in));
end

%% Add tide

if true

%Calculate tide
tidecon=zeros(length(tide.freq),4,length(out.lon));
for i3=1:size(tidecon,1)
    tidecon(i3,1,:)=tide.major_interp{i3}(out.lon,out.lat);
    tidecon(i3,3,:)=tide.phase_interp{i3}(out.lon,out.lat);
end
tidecon(:,3:4,:)=rad2deg(angle(tidecon(:,3:4,:)));
tt=struct('name',tide.name,'freq',tide.freq,'tidecon',zeros(length(tide.freq),4),'type','nodal');

if strcmpi(assim_type,'instant')
    for i1=1:length(out.ssh)
        tt.tidecon=squeeze(tidecon(:,:,i1));
        %add tide to output
        out.ssh(i1)=out.ssh(i1)+t_predic(out.t(i1),tt,'latitude',out.lat(i3));
    end
elseif strcmpi(assim_type,'daily_average')
    tModel=tLim(1)+[0:1:round(diff(tLim)*24)]/24; 
      
    for i1=1:length(out.ssh)
        %Get time
        tTmp=unique([tModel,out.t(i1)-.5,out.t(i1)+.5]); 
        tTmp=tTmp(tTmp>=max(tLim(1),out.t(i1)-.5)&tTmp<=min(tLim(2),out.t(i1)+.5)); 
        tt.tidecon=squeeze(tidecon(:,:,i1)); 
        %add tide to output
        sshTmp=t_predic(tTmp,tt,'latitude',out.lat(i3));
        out.ssh(i1)=out.ssh(i1)+trapz(tTmp,sshTmp);
    end
else
    error('unknown assim_type');
end

end

%% remove small cycles and demean

[cUni,c2uni,uni2c]=unique(out.cycle); 
for iCycle=1:length(cUni)
    in=uni2c==iCycle & ~isnan(out.ssh); 
    if sum(in)<cycle_length
        out.ssh(in)=NaN; 
    else
        out.ssh(in)=out.ssh(in)-mean(out.ssh(in)); 
    end
end
    

%% create output

%write output
in=~isnan(out.ssh); 
K=sum(in); 
Ktotal=length(ncread(outfile,'obs')); 
timeRef=tLim(1); 
if K>0
 ncwrite(outfile,'obs',out.ssh(in),Ktotal+1);
 ncwrite(outfile,'sig_d',out.sig(in),Ktotal+1);
 ncwrite(outfile,'type',out.cycle(in),Ktotal+1);
 ncwrite(outfile,'lon',out.lon(in),Ktotal+1);
 ncwrite(outfile,'lat',out.lat(in),Ktotal+1);
 ncwrite(outfile,'z',zeros(K,1),Ktotal+1);
 ncwrite(outfile,'time',(out.t(in)-timeRef)*24*3600,Ktotal+1);
end



end

%% read first and last time in file from its nam

function [tB,tE]=jason_tLim(fpath)

[fdir,fname,fext]=fileparts(fpath); 
line=regexp(fname,'(.*)_(.*)_(.*)_(.*)_(\d*)_(\d*)_(\d*)_(\d*)','tokens'); 
line=line{1};

tB=datenum([line{5},' ',line{6}],'yyyymmdd HHMMSS'); 
tE=datenum([line{7},' ',line{8}],'yyyymmdd HHMMSS'); 

end

function [tB,tE]=cryosat_tLim(fpath)

[fdir,fname,fext]=fileparts(fpath); 
line=regexp(fname,'c2p(\d+)c(\d+)_(\d*)_(\d*)','tokens'); 
line=line{1}; 

tB=datenum(line{3},'yyyymmdd'); 
tE=datenum(line{4},'yyyymmdd'); 

end

%% Median filter to remove outliers

function ssh=outlier_median(nsig,lon,lat,ssh)

for i1=1:length(ssh)
    d=deg2rad(hypot((lon-lon(i1))*cos(deg2rad(lat(i1))),lat-lat(i1)))*earthRadius;
    in1=d<=20e3 & d>0;
    if ~any(in1); ssh(i1)=NaN; continue; end
    
    median_ssh=nanmedian(ssh(in1));
    std_ssh=nanstd(ssh(in1));
    if abs(ssh(i1)-median_ssh)>nsig*std_ssh; ssh(i1)=NaN; end
    
end

end

%% Bayes filter to remove outliers

function [ssh,p]=outlier_bayes(sig,p_tresh,ssh)
    
    p_out=.02; 
    p_uni=1/(nanmax(ssh)-nanmin(ssh)); 
    ssh=ssh(:); sig=sig(:); 
    
    %Differences 
    dL=[NaN;diff(ssh)];
    dR=[dL(2:end);NaN]; 
    pL=pdf('norm',dL,zeros(size(dL)),[NaN;sqrt(sig(1:end-1).^2+sig(2:end).^2)]); 
    pR=pdf('norm',dR,zeros(size(dL)),[sqrt(sig(1:end-1).^2+sig(2:end).^2);NaN]); 

    %P(out|d1,d2)=P(d1,d2|out)P(out)/[P(d1,d2|out)P(out)+P(d1,d2|non-out)P(non-out)]
    p=nan(size(ssh)); 
    in=~isnan(dL) & ~isnan(dR); 
    p(in)=p_uni*p_uni*p_out./(p_uni*p_uni*p_out+pL(in).*pR(in)*(1-p_out)); 
    in=~isnan(dL) & isnan(dR); 
    p(in)=p_uni*p_out./(p_uni*p_out+pL(in)*(1-p_out)); 
    in=isnan(dL) & ~isnan(dR); 
    p(in)=p_uni*p_out./(p_uni*p_out+pR(in)*(1-p_out)); 
    
    %Set outliers to NaN
    in=p>p_tresh; 
     %display(sprintf('Removed %d outliers',sum(in))); 
    ssh(in)=NaN; 
    

end

