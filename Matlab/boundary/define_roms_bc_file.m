function [tmp]=define_roms_bc_file(filename,xi_rho,eta_rho,s_rho,...
    theta_s,theta_b,Tcline,Vtransform,Vstretching,grdname,timeunits,...
    north_yes,south_yes,west_yes,east_yes)

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

if east_yes
 nccreate(filename,'zeta_east','Dimensions',{'eta_rho',eta_rho,'time',Inf});
 nccreate(filename,'ubar_east','Dimensions',{'eta_u',eta_u,'time',Inf});
 nccreate(filename,'vbar_east','Dimensions',{'eta_v',eta_v,'time',Inf});
 nccreate(filename,'u_east','Dimensions',{'eta_u',eta_u,'s_rho',s_rho,'time',Inf});
 nccreate(filename,'v_east','Dimensions',{'eta_v',eta_v,'s_rho',s_rho,'time',Inf});
 nccreate(filename,'temp_east','Dimensions',{'eta_rho',eta_rho,'s_rho',s_rho,'time',Inf});
 nccreate(filename,'salt_east','Dimensions',{'eta_rho',eta_rho,'s_rho',s_rho,'time',Inf});
 
 ncwriteatt(filename,'zeta_east','time','ocean_time');
 ncwriteatt(filename,'ubar_east','time','ocean_time');
 ncwriteatt(filename,'vbar_east','time','ocean_time');
 ncwriteatt(filename,'u_east','time','ocean_time');
 ncwriteatt(filename,'v_east','time','ocean_time');
 ncwriteatt(filename,'temp_east','time','ocean_time');
 ncwriteatt(filename,'salt_east','time','ocean_time');


end

if west_yes
 nccreate(filename,'zeta_west','Dimensions',{'eta_rho',eta_rho,'time',Inf});
 nccreate(filename,'ubar_west','Dimensions',{'eta_u',eta_u,'time',Inf});
 nccreate(filename,'vbar_west','Dimensions',{'eta_v',eta_v,'time',Inf});
 nccreate(filename,'u_west','Dimensions',{'eta_u',eta_u,'s_rho',s_rho,'time',Inf});
 nccreate(filename,'v_west','Dimensions',{'eta_v',eta_v,'s_rho',s_rho,'time',Inf});
 nccreate(filename,'temp_west','Dimensions',{'eta_rho',eta_rho,'s_rho',s_rho,'time',Inf});
 nccreate(filename,'salt_west','Dimensions',{'eta_rho',eta_rho,'s_rho',s_rho,'time',Inf});

 ncwriteatt(filename,'zeta_west','time','ocean_time');
 ncwriteatt(filename,'ubar_west','time','ocean_time');
 ncwriteatt(filename,'vbar_west','time','ocean_time');
 ncwriteatt(filename,'u_west','time','ocean_time');
 ncwriteatt(filename,'v_west','time','ocean_time');
 ncwriteatt(filename,'temp_west','time','ocean_time');
 ncwriteatt(filename,'salt_west','time','ocean_time');

end




if north_yes
 nccreate(filename,'zeta_north','Dimensions',{'xi_rho',xi_rho,'time',Inf});
 nccreate(filename,'ubar_north','Dimensions',{'xi_u',xi_u,'time',Inf});
 nccreate(filename,'vbar_north','Dimensions',{'xi_v',xi_v,'time',Inf});
 nccreate(filename,'u_north','Dimensions',{'xi_u',xi_u,'s_rho',s_rho,'time',Inf});
 nccreate(filename,'v_north','Dimensions',{'xi_v',xi_v,'s_rho',s_rho,'time',Inf});
 nccreate(filename,'temp_north','Dimensions',{'xi_rho',xi_rho,'s_rho',s_rho,'time',Inf});
 nccreate(filename,'salt_north','Dimensions',{'xi_rho',xi_rho,'s_rho',s_rho,'time',Inf});

 ncwriteatt(filename,'zeta_north','time','ocean_time');
 ncwriteatt(filename,'ubar_north','time','ocean_time');
 ncwriteatt(filename,'vbar_north','time','ocean_time');
 ncwriteatt(filename,'u_north','time','ocean_time');
 ncwriteatt(filename,'v_north','time','ocean_time');
 ncwriteatt(filename,'temp_north','time','ocean_time');
 ncwriteatt(filename,'salt_north','time','ocean_time');
 
end

if south_yes
 nccreate(filename,'zeta_south','Dimensions',{'xi_rho',xi_rho,'time',Inf});
 nccreate(filename,'ubar_south','Dimensions',{'xi_u',xi_u,'time',Inf});
 nccreate(filename,'vbar_south','Dimensions',{'xi_v',xi_v,'time',Inf});
 nccreate(filename,'u_south','Dimensions',{'xi_u',xi_u,'s_rho',s_rho,'time',Inf});
 nccreate(filename,'v_south','Dimensions',{'xi_v',xi_v,'s_rho',s_rho,'time',Inf});
 nccreate(filename,'temp_south','Dimensions',{'xi_rho',xi_rho,'s_rho',s_rho,'time',Inf});
 nccreate(filename,'salt_south','Dimensions',{'xi_rho',xi_rho,'s_rho',s_rho,'time',Inf});

 ncwriteatt(filename,'zeta_south','time','ocean_time');
 ncwriteatt(filename,'ubar_south','time','ocean_time');
 ncwriteatt(filename,'vbar_south','time','ocean_time');
 ncwriteatt(filename,'u_south','time','ocean_time');
 ncwriteatt(filename,'v_south','time','ocean_time');
 ncwriteatt(filename,'temp_south','time','ocean_time');
 ncwriteatt(filename,'salt_south','time','ocean_time');

end

ncwriteatt(filename,'ocean_time','units',timeunits);


ncwriteatt(filename,'/','type','ROMS BOUNDARY file');
ncwriteatt(filename,'/','grd_name',grdname);

ncwrite(filename,'theta_s',theta_s);
ncwrite(filename,'theta_b',theta_b);
ncwrite(filename,'Tcline',Tcline);
ncwrite(filename,'Vtransform',Vtransform);
ncwrite(filename,'Vstretching',Vstretching);

tmp=0;
