function [tmp]=define_roms_cloud_file(filename,xi_rho,eta_rho,grdname,timeunits)

%xi_u=xi_rho-1;
%eta_u=eta_rho;
%xi_v=xi_rho;
%eta_v=eta_rho-1;

if exist(filename,'file')==2
   eval(['!rm ' filename]);
end

nccreate(filename,'cloud_time','Dimensions',{'time',Inf},'format','classic');
nccreate(filename,'cloud','Dimensions',...
  {'xi_rho',xi_rho,'eta_rho',eta_rho,'time',Inf});

ncwriteatt(filename,'cloud_time','units',timeunits);
ncwriteatt(filename,'cloud_time','long_name','cloud fraction time');

ncwriteatt(filename,'cloud','long_name','cloud fraction (0 to 1)');
ncwriteatt(filename,'cloud','units','nondimensional');
ncwriteatt(filename,'cloud','time','cloud_time');

ncwriteatt(filename,'/','type','ROMS FORCING file');
ncwriteatt(filename,'/','grd_name',grdname);

tmp=0;