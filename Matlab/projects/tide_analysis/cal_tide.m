%% calculate tides
clear all; clc; 
addpath('/home/server/student/homes/ipasmans/Matlab/roms/t_tide/'); 
addpath('/home/server/student/homes/ipasmans/Matlab/figures/'); 

roms_dir={'/home/aruba/vol2/ipasmans/Exp/Exp26_free/'}; 
grd_file='/home/aruba/vol2/ipasmans/Exp/Prm/grd_ow2km_r13lin_mix.nc'; 
out_file='tideUV_Exp26fstride1'; 
t=[datenum('2011-04-01'):1/24:datenum('2011-10-01')];
dateRef=datenum('2005-01-01'); 


%% read grid

grd.lon=ncread(grd_file,'lon_rho'); 
grd.lat=ncread(grd_file,'lat_rho'); 
grd.mask=ncread(grd_file,'mask_rho'); 
grd.sig=ncread('/home/aruba/vol2/ipasmans/Exp/Prm/per_sig.nc','temp'); 

%% make output

tide.name=['K2';'S2';'M2';'N2';'K1';'P1';'O1';'Q1']; 

out.lon=grd.lon; 
out.lat=grd.lat; 
out.val=nan(size(grd.lon,1),size(grd.lon,2),size(tide.name,1)); 
out.dval=nan(size(grd.lon,1),size(grd.lon,2),size(tide.name,1)); 
out.freq=[]; 
out.name=[]; 
out.t=t; 

%% calculate tides

stride=1;
zeta=nan(ceil(size(grd.lon,1)/stride),ceil(size(grd.lon,2)/stride),length(t));

for iDir=1:length(roms_dir)
    con=dir(roms_dir{iDir});
    for iFile=1:length(con)
        
        if strncmpi(con(iFile).name,'ocean_his_',10)
            filename=fullfile(roms_dir{iDir},con(iFile).name); 
            display(filename); 
            
            tFile=ncread(filename,'ocean_time');
            tFile=tFile/24/3600+dateRef;
            [int,it]=ismember(tFile,t);
            if ~any(int); continue; end
            
            it1=find(int,1,'first');
            it2=find(int,1,'last');
            nt=it2-it1+1;
            zeta1=ncread(filename,'zeta',[1,1,it1],[Inf,Inf,nt],[stride,stride,1]);
           
            
            
            tFile=tFile(it1:it2);
            for k=1:length(tFile)
                it3=find(t==tFile(k),1,'first'); 
                if isempty(it3); continue; end
                zeta(:,:,it3)=zeta1; 
            end
            
        end
    end
end
   

for i2=[1:stride:size(grd.lon,2)]
    for i1=[1:stride:size(grd.lon,1)]
        if ~grd.mask(i1,i2); continue; end
        %if grd.sig(i1,i2)<max(grd.sig(:)); continue; end
        
        display(sprintf('Calculating point (%d,%d)',i1,i2)); 
        
        zeta1=squeeze(zeta( (i1-1)/stride+1,(i2-1)/stride+1,:)); 

        tide=t_tide(zeta1,'interval',1,'start time',t(1),'latitude',grd.lat(i1,i2),'secular','mean',...
            'Rayleigh',tide.name,'output','none'); 
        
        out.name=tide.name; 
        out.freq=tide.freq; 
        out.val(i1,i2,:)=tide.tidecon(:,1).*exp(1i*deg2rad(tide.tidecon(:,3))); 
        out.dval(i1,i2,:)=tide.tidecon(:,2).*exp(1i*deg2rad(tide.tidecon(:,4))); 
        
    end
end

%% interpolant

% for i3=1:size(out.val,3)
%     
% for i2=1:size(grd.lat,2)
%    inan=~isnan(out.val(:,i2,i3)); 
%    if sum(inan)>=2
%        out.val(:,i2,i3)=interp1(grd.lon(inan,i2),out.val(inan,i2,i3),grd.lon(:,i2)); 
%        out.dval(:,i2,i3)=interp1(grd.lon(inan,i2),out.dval(inan,i2,i3),grd.lon(:,i2));
%    end
% end
% 
% for i1=1:size(grd.lon,1)
%     inan=~isnan(out.val(i1,:,i3)); 
%     if sum(inan)>=2
%         out.val(i1,:,i3)=interp1(grd.lat(i1,inan),out.val(i1,inan,i3),grd.lat(i1,:)); 
%         out.dval(i1,:,i3)=interp1(grd.lat(i1,inan),out.dval(i1,inan,i3),grd.lat(i1,:)); 
%     end
% end
% 
%     
% end



%% 

save(out_file,'-struct','out'); 
