function [tmp]=define_roms_pair_file(filename,xi_rho,eta_rho,grdname,timeunits)

%xi_u=xi_rho-1;
%eta_u=eta_rho;
%xi_v=xi_rho;
%eta_v=eta_rho-1;

if exist(filename,'file')==2
   eval(['!rm ' filename]);
end

nccreate(filename,'pair_time','Dimensions',{'time',Inf},'format','classic');
nccreate(filename,'Pair','Dimensions',...
  {'xi_rho',xi_rho,'eta_rho',eta_rho,'time',Inf});

ncwriteatt(filename,'pair_time','units',timeunits);
ncwriteatt(filename,'pair_time','long_name','surface air pressure time');

ncwriteatt(filename,'Pair','long_name','surface air pressure');
ncwriteatt(filename,'Pair','units','milibar');
ncwriteatt(filename,'Pair','time','pair_time');

ncwriteatt(filename,'/','type','ROMS FORCING file');
ncwriteatt(filename,'/','grd_name',grdname);

tmp=0;