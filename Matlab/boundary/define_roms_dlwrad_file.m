function [tmp]=define_roms_dlwrad_file(filename,xi_rho,eta_rho,grdname,timeunits)

%xi_u=xi_rho-1;
%eta_u=eta_rho;
%xi_v=xi_rho;
%eta_v=eta_rho-1;

if exist(filename,'file')==2
   eval(['!rm ' filename]);
end

nccreate(filename,'lrf_time','Dimensions',{'time',Inf},'format','classic');
nccreate(filename,'lwrad_down','Dimensions',...
  {'xi_rho',xi_rho,'eta_rho',eta_rho,'time',Inf});

ncwriteatt(filename,'lrf_time','units',timeunits);
ncwriteatt(filename,'lrf_time','long_name','longwave radiation time');

ncwriteatt(filename,'lwrad_down','long_name','downward longwave radiation');
ncwriteatt(filename,'lwrad_down','units','W m-2');
ncwriteatt(filename,'lwrad_down','positive value','downward flux, heating');
ncwriteatt(filename,'lwrad_down','time','lrf_time');

ncwriteatt(filename,'/','type','ROMS FORCING file');
ncwriteatt(filename,'/','grd_name',grdname);

tmp=0;