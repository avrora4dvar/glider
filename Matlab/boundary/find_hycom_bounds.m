function [i1,i2,j1,j2,xb,yb,XX,YY]=find_hycom_bounds(lon,lat,X,Y,Lon,Lat)

% Find bounding indices for X and Y:
% USAGE: [i1,i2,j1,j2,xb,yb,XX,YY]=find_hycom_bounds(lons,lats,X,Y,Lon,Lat)
% INPUTS: lon = 2D, roms subgrid longitudes (-180 180)
%         lat = 2D, roms subgrid latitudes  bounds (-90  90)
%         X = vector of HYCOM x-coordinates
%         Y = vector of HYCOM y-coordinates
%         Lon = a 2D array of HYCOM longitudes (-180:180)
%         Lat = a 2D array of HYCOM latitudes (-90 90)
%
% OUTPUTS: i1,i2,j1,j2: limit indices of the HYCOM box surrounding the
%                       roms (lon,lat) domain
%          xb, yb: coordinates of the roms grid in terms of HYCOM indices
%          XX, YY: index coordinates of the HYCOM grid produced by meshgrid
%               ( note: ii=i1:i2; jj=j1:j2, [YY,XX]=meshgrid(Y(jj),X(ii)) )

lons(1)=min(min(lon));
lons(2)=max(max(lon));
lats(1)=min(min(lat));
lats(2)=max(max(lat));

[nx,ny]=size(Lon);
n=nx*ny;
[YY,XX]=meshgrid(Y,X);
Lon1D=reshape(Lon,[n 1]);
Lat1D=reshape(Lat,[n 1]);
X1D=reshape(XX,[n 1]);
Y1D=reshape(YY,[n 1]);

[yyl,xxl]=meshgrid(lats,lons);
xl=reshape(xxl,[4 1]);
yl=reshape(yyl,[4 1]);

Xk=zeros(4,1);
Yk=zeros(4,1);
for k=1:4
 [tmp,i1]=min((Lon1D-xl(k)).^2+(Lat1D-yl(k)).^2);
 Xk(k)=X1D(i1);
 Yk(k)=Y1D(i1);
end

i1=min(Xk)-2;
i2=max(Xk)+2;
j1=min(Yk)-2;
j2=max(Yk)+2;

ii=[i1:i2];
jj=[j1:j2];
nx=length(ii);
ny=length(jj);
[YY,XX]=meshgrid(Y(jj),X(ii));

% roms coordinates in terms of the HYCOM X and Y:
dlon=1/25;
dlat=1/25;
[lat2,lon2]=meshgrid(...
    floor(min(min(Lat(ii,jj)))):dlat:ceil(max(max(Lat(ii,jj)))),...
    floor(min(min(Lon(ii,jj)))):dlon:ceil(max(max(Lon(ii,jj)))));
                 
Xgridded=griddata(reshape(Lat(ii,jj),[nx*ny 1]),...
                  reshape(Lon(ii,jj),[nx*ny 1]),...
                  reshape(XX,[nx*ny 1]),...
                  lat2,lon2);
Ygridded=griddata(reshape(Lat(ii,jj),[nx*ny 1]),...
                  reshape(Lon(ii,jj),[nx*ny 1]),...
                  reshape(YY,[nx*ny 1]),...
                  lat2,lon2);
xb=interp2(lat2,lon2,Xgridded,lat,lon);
yb=interp2(lat2,lon2,Ygridded,lat,lon);
