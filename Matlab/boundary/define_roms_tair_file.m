function [tmp]=define_roms_tair_file(filename,xi_rho,eta_rho,grdname,timeunits)

%xi_u=xi_rho-1;
%eta_u=eta_rho;
%xi_v=xi_rho;
%eta_v=eta_rho-1;

if exist(filename,'file')==2
   eval(['!rm ' filename]);
end

nccreate(filename,'tair_time','Dimensions',{'time',Inf},'format','classic');
nccreate(filename,'Tair','Dimensions',...
  {'xi_rho',xi_rho,'eta_rho',eta_rho,'time',Inf});

ncwriteatt(filename,'tair_time','units',timeunits);
ncwriteatt(filename,'tair_time','long_name','surface air temperature time');

ncwriteatt(filename,'Tair','long_name','surface air temperature');
ncwriteatt(filename,'Tair','units','Celsius');
ncwriteatt(filename,'Tair','time','tair_time');

ncwriteatt(filename,'/','type','ROMS FORCING file');
ncwriteatt(filename,'/','grd_name',grdname);

tmp=0;
