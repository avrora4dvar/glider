%% Plot covariances
addpath('../io'); 
addpath('../../projects/obslist'); 

%% Read grid

grdFile='/home/aruba/vol2/ipasmans/Exp/Prm/grd_ow2km_r13lin_mix.nc';
grd.lon=ncread(grdFile,'lon_rho'); 
grd.lat=ncread(grdFile,'lat_rho'); 
grd.h=ncread(grdFile,'h'); 
grd.zr=ncread(grdFile,'z0_r'); 
grd.zw=ncread(grdFile,'z0_w'); 

%% Read obslist

obs=io_glider_osu('~/Obs/glider/osu','dir'); 
obs.t=repmat(obs.t,[1,size(obs.temp,2)]); 
obs.profile.z=[2:4:248]; 
for iz=1:length(obs.profile.z)
    in=obs.depth>=obs.profile.z(iz)-2 & obs.depth<obs.profile.z(iz)+2; 
    in=in & obs.t>=datenum('2011-07-21') & obs.t<=datenum('2011-08-14'); 
    in=in(:); 
    obs.profile.temp(iz)=nanmean(obs.temp(in)); 
    obs.profile.salt(iz)=nanmean(obs.salt(in));
end
temp1=interp1([min(obs.profile.temp),max(obs.profile.temp)],[0,1],obs.profile.temp); 

%%

obsDir='/home/aruba/vol2/ipasmans/Exp/Obs/Exp37'; 
obs.profile.sig_temp=zeros(size(obs.profile.z)); 
obs.profile.sig_N=zeros(size(obs.profile.z)); 

for iWin=[2392:3:2410]
   fname=fullfile(obsDir,sprintf('obslist%-.4d.nc',iWin)); 
   obs1=read_obslist(fname); 
   for iz=1:length(obs.profile.z)
    in=obs1.type==6 & obs1.z>=obs.profile.z(iz)-2 & obs1.z<obs.profile.z(iz)+2;
    if ~any(in); continue; end
        obs.profile.sig_temp(iz)=obs.profile.sig_temp(iz)+sum(obs1.sig(in).^2); 
        obs.profile.sig_N(iz)=obs.profile.sig_N(iz)+sum(in); 
   end
end
obs.profile.sig_temp=obs.profile.sig_temp./obs.profile.sig_N; 

%% Calculate 

z1=squeeze(grd.zr(203,132,:)); 
c1=.9^2*exp(z1/100).*exp(-.5*(z1-0).^2/15^2); 
c2=.9^2*exp(z1/100).*exp(-.5*(z1+100).^2/15^2); 

%% Save

profile=obs.profile; 
save('~/Data/article_glider/cov_z.mat','profile','c1','c2','z1'); 

%% Plot

clf; hold on;
plot(temp1,-obs.profile.z,'k-'); 
plot(c1,z1,'b-'); 
plot(c2,z1,'g-'); 
plot(obs.profile.sig_temp,-obs.profile.z,'r--'); 
