%% calculate tides
clear all; clc; 
addpath('/home/server/student/homes/ipasmans/Matlab/roms/t_tide/'); 
addpath('/home/server/student/homes/ipasmans/Matlab/figures/'); 

roms_dir={'/home/aruba/vol2/ipasmans/Exp/Exp26_free/'}; 
grd_file='/home/aruba/vol2/ipasmans/Exp/Prm/grd_ow2km_r13lin_mix.nc'; 
out_file='tideUV40_Exp26fstride1'; 
t=[datenum('2011-04-01'):1/24:datenum('2011-10-01')];
dateRef=datenum('2005-01-01'); 
stride=1;


%% read grid

grd.lon=ncread(grd_file,'lon_rho'); 
grd.lat=ncread(grd_file,'lat_rho'); 
grd.mask=ncread(grd_file,'mask_rho'); 
grd.sig=ncread('/home/aruba/vol2/ipasmans/Exp/Prm/per_sig_bal.nc','temp'); 

%% make output

tide.name=['K2';'S2';'M2';'N2';'K1';'P1';'O1';'Q1']; 

out.lon=grd.lon; 
out.lat=grd.lat; 
out.major=nan(size(grd.lon,1),size(grd.lon,2),size(tide.name,1)); 
out.minor=out.major; out.dir=out.major; out.phase=out.major; 
out.dmajor=nan(size(grd.lon,1),size(grd.lon,2),size(tide.name,1)); 
out.dminor=out.dmajor; out.ddir=out.dmajor; out.dphase=out.dmajor;
out.freq=[]; 
out.name=[]; 
out.t=t; 

s=size(grd.lon); 
countu=min(length([1:stride:size(grd.lon,1)-1]),length([2:stride:size(grd.lon,1)-1])); 
countv=min(length([1:stride:size(grd.lon,2)-1]),length([2:stride:size(grd.lon,2)-1])); 

%% calculate tides


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
            
            
%             uv1=.5*ncread(filename,'ubar',[1,2,it1],[countu,countv,nt],[stride,stride,1])...
%                 +.5*ncread(filename,'ubar',[2,2,it1],[countu,countv,nt],[stride,stride,1])...
%                 +.5i*ncread(filename,'vbar',[2,1,it1],[countu,countv,nt],[stride,stride,1])...
%                 +.5i*ncread(filename,'vbar',[2,2,it1],[countu,countv,nt],[stride,stride,1]);   

            uv1=.5*ncread(filename,'u',[1,2,40,it1],[countu,countv,1,nt],[stride,stride,1,1])...
                +.5*ncread(filename,'u',[2,2,40,it1],[countu,countv,1,nt],[stride,stride,1,1])...
                +.5i*ncread(filename,'v',[2,1,40,it1],[countu,countv,1,nt],[stride,stride,1,1])...
                +.5i*ncread(filename,'v',[2,2,40,it1],[countu,countv,1,nt],[stride,stride,1,1]);     
        
            
            tFile=tFile(it1:it2);
            for k=1:length(tFile)
                it3=find(t==tFile(k),1,'first'); 
                if isempty(it3); continue; end
                uv(:,:,it3)=uv1(:,:,k); 
            end
            
        end
    end
end
   

for i2=[2:stride:size(grd.lon,2)-1]
    for i1=[2:stride:size(grd.lon,1)-1]
        if ~grd.mask(i1,i2); continue; end
        %if grd.sig(i1,i2)<max(grd.sig(:)); continue; end
        
        display(sprintf('Calculating point (%d,%d)',i1,i2)); 
        
        uv1=squeeze(uv( (i1-2)/stride+1,(i2-2)/stride+1,:)); 
        if ~any(~isnan(uv1)); continue; end

        u1=real(uv1(:)); v1=imag(uv1(:)); 
        u1=u1-filter_lanczos(u1,1,1/40,60,1); 
        v1=v1-filter_lanczos(v1,1,1/40,60,1); 
        
        tide=t_tide(u1+1i*v1,'interval',1,'start time',t(1),'latitude',grd.lat(i1,i2),'secular','mean',...
            'Rayleigh',tide.name,'output','none'); 
        
        out.name=tide.name; 
        out.freq=tide.freq; 
        %Values
        out.major(i1,i2,:)=tide.tidecon(:,1);
        out.minor(i1,i2,:)=tide.tidecon(:,3);
        out.dir(i1,i2,:)=tide.tidecon(:,5); 
        out.phase(i1,i2,:)=tide.tidecon(:,7); 
        %Errors
        out.dmajor(i1,i2,:)=tide.tidecon(:,2); 
        out.dminor(i1,i2,:)=tide.tidecon(:,4); 
        out.ddir(i1,i2,:)=tide.tidecon(:,6);
        out.dphase(i1,i2,:)=tide.tidecon(:,8); 
        
   
        
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
