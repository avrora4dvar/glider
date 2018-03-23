% replaces time to cover year 2014 in existing CR file - quick fix 
% before making near-real time upgrade script
fname_in='../Prm/ow2_cr4_offset0.nc';
fname_out='../Prm/ow2km_CR_N40.nc';
t0=ncread(fname_in,'river_time');
dates=datestr(datenum(t0+datenum(2005,1,1)));
% 
shift_days=datenum(2012,1,1,12,0,0)-datenum(dates(1,:));
t=t0+shift_days;
ncwrite(fname_out,'river_time',t);
 
