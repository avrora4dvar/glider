function [tmp]=define_roms_wind_file(filename,xi_rho,eta_rho,grdname,timeunits)

%xi_u=xi_rho-1;
%eta_u=eta_rho;
%xi_v=xi_rho;
%eta_v=eta_rho-1;

if exist(filename,'file')==2
   eval(['!rm ' filename]);
end

nccreate(filename,'wind_time','Dimensions',{'time',Inf},'format','classic');
nccreate(filename,'Uwind','Dimensions',...
  {'xi_rho',xi_rho,'eta_rho',eta_rho,'time',Inf});
nccreate(filename,'Vwind','Dimensions',...
  {'xi_rho',xi_rho,'eta_rho',eta_rho,'time',Inf});

ncwriteatt(filename,'wind_time','units',timeunits);
ncwriteatt(filename,'wind_time','long_name','surface wind time');

ncwriteatt(filename,'Uwind','long_name','surface u-wind component');
ncwriteatt(filename,'Uwind','units','m s-1');
ncwriteatt(filename,'Uwind','time','wind_time');

ncwriteatt(filename,'Vwind','long_name','surface v-wind component');
ncwriteatt(filename,'Vwind','units','m s-1');
ncwriteatt(filename,'Vwind','time','wind_time');

ncwriteatt(filename,'/','type','ROMS FORCING file');
ncwriteatt(filename,'/','grd_name',grdname);

tmp=0;
