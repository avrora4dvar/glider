function [tmp]=define_roms_ic_file(filename,xi_rho,eta_rho,s_rho,...
    theta_s,theta_b,Tcline,Vtransform,Vstretching,grdname,timeunits)

xi_u=xi_rho-1;
eta_u=eta_rho;
xi_v=xi_rho;
eta_v=eta_rho-1;

if exist(filename,'file')==2
   eval(['!rm ' filename]);
end

nccreate(filename,'theta_s','format','classic');
nccreate(filename,'theta_b');
nccreate(filename,'Tcline');
nccreate(filename,'Vtransform');
nccreate(filename,'Vstretching');
nccreate(filename,'ocean_time','Dimensions',{'time',Inf});

nccreate(filename,'zeta','Dimensions',{'xi_rho',xi_rho,'eta_rho',eta_rho,'time',Inf});
nccreate(filename,'ubar','Dimensions',{'xi_u',xi_u,'eta_u',eta_u,'time',Inf});
nccreate(filename,'vbar','Dimensions',{'xi_v',xi_v,'eta_v',eta_v,'time',Inf});
nccreate(filename,'u','Dimensions',{'xi_u',xi_u,'eta_u',eta_u,'s_rho',s_rho,'time',Inf});
nccreate(filename,'v','Dimensions',{'xi_v',xi_v,'eta_v',eta_v,'s_rho',s_rho,'time',Inf});
nccreate(filename,'temp','Dimensions',{'xi_rho',xi_rho,'eta_rho',eta_rho,'s_rho',s_rho,'time',Inf});
nccreate(filename,'salt','Dimensions',{'xi_rho',xi_rho,'eta_rho',eta_rho,'s_rho',s_rho,'time',Inf});

ncwriteatt(filename,'ocean_time','units',timeunits);
ncwriteatt(filename,'zeta','time','ocean_time');
ncwriteatt(filename,'ubar','time','ocean_time');
ncwriteatt(filename,'vbar','time','ocean_time');
ncwriteatt(filename,'u','time','ocean_time');
ncwriteatt(filename,'v','time','ocean_time');
ncwriteatt(filename,'temp','time','ocean_time');
ncwriteatt(filename,'salt','time','ocean_time');

ncwriteatt(filename,'/','type','ROMS INITIAL file');
ncwriteatt(filename,'/','grd_name',grdname);

ncwrite(filename,'theta_s',theta_s);
ncwrite(filename,'theta_b',theta_b);
ncwrite(filename,'Tcline',Tcline);
ncwrite(filename,'Vtransform',Vtransform);
ncwrite(filename,'Vstretching',Vstretching);

tmp=0;