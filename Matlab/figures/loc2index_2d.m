function [ix,iy,w] = loc2index_2d(grd,loc)
%LOC2INDEX_LONLAT [ix,iy,w]=loc2index_2d
%Get grid indices and weighting factors for bilinear
%interpolation. Here loc.lon,loc.lat are 1xn matrices with longitude and
%latitude and ix,iy,w are 4xn matrices with indices and weights. 

if ~isfield(grd,'mask')
    grd.mask=ones(size(grd.lon)); 
end

%initiate output
npoints=length(loc.lon); 
ix=nan(4,npoints); 
iy=nan(4,npoints); 
w=zeros(4,npoints); 

%Lon between 0-360
grd.lon=mod(grd.lon,360); 
loc.lon=mod(loc.lon,360); 

%Interpolate indices
ix1=interp1(reshape(grd.lon(:,1),[],1),[1:length(grd.lon(:,1))],loc.lon); 
iy1=interp1(reshape(grd.lat(1,:),[],1),[1:length(grd.lat(1,:))],loc.lat); 

%Assign index ix
ix(1,:)=floor(ix1); 
ix(2,:)=floor(ix1); 
ix(3,:)=ceil(ix1); 
ix(4,:)=ceil(ix1);

%Assign weight w
w(1,:)=1-abs(ix1-ix(1,:)); 
w(2,:)=1-abs(ix1-ix(2,:)); 
w(3,:)=1-abs(ix1-ix(3,:)); 
w(4,:)=1-abs(ix1-ix(4,:)); 

%Assign index iy
iy(1,:)=floor(iy1); 
iy(2,:)=ceil(iy1); 
iy(3,:)=floor(iy1); 
iy(4,:)=ceil(iy1); 

%Assign weight
w(1,:)=w(1,:).*(1-abs(iy1-iy(1,:))); 
w(2,:)=w(2,:).*(1-abs(iy1-iy(2,:))); 
w(3,:)=w(3,:).*(1-abs(iy1-iy(3,:))); 
w(4,:)=w(4,:).*(1-abs(iy1-iy(4,:))); 

%Mask out
w(isnan(w))=0;
w(isnan(ix)|isnan(iy))=0; 
for i2=1:size(w,2)
    for i1=1:size(w,1)
        if w(i1,i2)==0; continue; end
        w(i1,i2)=w(i1,i2)*grd.mask(ix(i1,i2),iy(i1,i2)); 
    end
end

%Normalize w
sumw=sum(w,1);
in=sumw>0; 
w(:,in)=bsxfun(@times,w(:,in),1./sumw(in)); 
           

end

