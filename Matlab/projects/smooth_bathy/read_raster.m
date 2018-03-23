%% read Puget Sound

fname='~/Data/bathymetry/psdem_2005/psdem_wgs.asc'; 
fid=fopen(fname); 
line=fgetl(fid); ncol=sscanf(line,'ncols %d'); 
line=fgetl(fid); nrows=sscanf(line,'nrows %d'); 
line=fgetl(fid); xll=sscanf(line,'xllcorner %f');
line=fgetl(fid); yll=sscanf(line,'yllcorner %f'); 
line=fgetl(fid); cellsize=sscanf(line,'cellsize %f'); 
line=fgetl(fid); nanVal=sscanf(line,'NODATA_value %f'); 
fclose(fid); 

x=xll+[0:ncol-1]*cellsize; 
y=yll+[nrows-1:-1:0]*cellsize; 
[x,y]=meshgrid(x,y); 

data=dlmread(fname,' ',6,0);
data(data==nanVal)=NaN; 
data=data(1:nrows,1:ncol); 

lon=x'; lat=y'; h=-double(data)'*.3048; 
meta='http://www.ocean.washington.edu/data/pugetsound/psdem2005.html';
save('~/Matlab/projects/smooth_bathy/puget_sound.mat','lon','lat','h','meta'); 

%% Read columbia River


fname='~/Data/bathymetry/P260_B90/Columbia3.asc'; 
fid=fopen(fname); 
line=fgetl(fid); ncol=sscanf(line,'ncols %d'); 
line=fgetl(fid); nrows=sscanf(line,'nrows %d'); 
line=fgetl(fid); xll=sscanf(line,'xllcorner %f');
line=fgetl(fid); yll=sscanf(line,'yllcorner %f'); 
line=fgetl(fid); cellsize=sscanf(line,'cellsize %f'); 
line=fgetl(fid); nanVal=sscanf(line,'NODATA_value %f'); 
fclose(fid); 

x=xll+[0:ncol-1]*cellsize; 
y=yll+[nrows-1:-1:0]*cellsize; 
[x,y]=meshgrid(x,y); 

data=dlmread(fname,' ',6,0);
data(data==nanVal)=NaN; 
data=data(1:nrows,1:ncol); 


lon=x'; lat=y'; h=-double(data)'+1.4; 
meta='http://estuarinebathymetry.noaa.gov/bathy_htmls/P260.html'; 
save('~/Matlab/projects/smooth_bathy/columbia_river3.mat','lon','lat','h','meta'); 