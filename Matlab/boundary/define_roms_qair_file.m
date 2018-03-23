function [tmp]=define_roms_qair_file(filename,xi_rho,eta_rho,grdname,timeunits)

%xi_u=xi_rho-1;
%eta_u=eta_rho;
%xi_v=xi_rho;
%eta_v=eta_rho-1;

if exist(filename,'file')==2
   eval(['!rm ' filename]);
end

nccreate(filename,'qair_time','Dimensions',{'time',Inf},'format','classic');
nccreate(filename,'Qair','Dimensions',...
  {'xi_rho',xi_rho,'eta_rho',eta_rho,'time',Inf});

ncwriteatt(filename,'qair_time','units',timeunits);
ncwriteatt(filename,'qair_time','long_name','surface relative humidity time');

ncwriteatt(filename,'Qair','long_name','surface relative humidity');
ncwriteatt(filename,'Qair','units','milibar');
ncwriteatt(filename,'Qair','time','qair_time');

ncwriteatt(filename,'/','type','ROMS FORCING file');
ncwriteatt(filename,'/','grd_name',grdname);

tmp=0;