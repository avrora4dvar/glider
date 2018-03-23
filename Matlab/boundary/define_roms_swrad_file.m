function [tmp]=define_roms_swrad_file(filename,xi_rho,eta_rho,grdname,timeunits)

%xi_u=xi_rho-1;
%eta_u=eta_rho;
%xi_v=xi_rho;
%eta_v=eta_rho-1;

if exist(filename,'file')==2
   eval(['!rm ' filename]);
end

nccreate(filename,'srf_time','Dimensions',{'time',Inf},'format','classic');
nccreate(filename,'swrad','Dimensions',...
  {'xi_rho',xi_rho,'eta_rho',eta_rho,'time',Inf});

ncwriteatt(filename,'srf_time','units',timeunits);
ncwriteatt(filename,'srf_time','long_name','solar shortwave radiation time');

ncwriteatt(filename,'swrad','long_name','solar shortwave radiation');
ncwriteatt(filename,'swrad','units','W m-2');
ncwriteatt(filename,'swrad','time','srf_time');
ncwriteatt(filename,'swrad','positive value','downward flux, heating');
ncwriteatt(filename,'swrad','negative value','upward flux, cooling');

ncwriteatt(filename,'/','type','ROMS FORCING file');
ncwriteatt(filename,'/','grd_name',grdname);

tmp=0;