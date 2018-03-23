function [ out ] = get_roms_2d(roms_dir,field_name,grd,loc,varargin)
%GET_ROMS_2D
% out=get_roms_2d(roms_dir,field_name,grd,loc,layer) reads field
% field_name from history files in roms_dir at locations in loc
% Here grd is struct with fields lon,lat,mask containing the grid
% associated with field_name and loc a struct with fields lon,lat,t
% containing the longitudes,latitudes,times (in datenum) where
% field_name is requested. layer is the s-layer from which output is
% taken (optional). By default surface layer is chosen. 

if length(varargin)==1
    layer=varargin{1}; nlayer=1; 
else
    layer=[]; nlayer=[]; 
end
npoints=length(loc.lon);

%get grid indices
[ix,iy,w]=loc2index_2d(grd,loc);
iBnd.ix(1)=nanmin(ix(:)); iBnd.ix(2)=nanmax(ix(:)); 
iBnd.ix(3)=iBnd.ix(2)-iBnd.ix(1)+1; 
iBnd.iy(1)=nanmin(iy(:)); iBnd.iy(2)=nanmax(iy(:)); 
iBnd.iy(3)=iBnd.iy(2)-iBnd.iy(1)+1;
ix=ix-iBnd.ix(1)+1; 
iy=iy-iBnd.iy(1)+1; 
ix(isnan(ix))=1; iy(isnan(iy))=1; 

%get roms indices
[hash,t_roms]=get_hash_roms(roms_dir,[min(loc.t),max(loc.t)]);
info=ncinfo(hash(1).filename,field_name); 
if length(info.Size)==4 && isempty(layer)
    layer=info.Size(3); nlayer=1; %set surface layer as default
end

%Time index
[it,wt]=loc2index_1d(t_roms,loc.t);

out=struct('lon',loc.lon,'lat',loc.lat,'t',loc.t,'field',field_name,'val',zeros(size(loc.t))); 
w_sum=zeros(size(loc.t)); 
for i0=1:length(hash)
    %Read from netcdf file
    ncstart=[iBnd.ix(1),iBnd.iy(1),layer,hash(i0).index(1)];
    nccount=[iBnd.ix(3),iBnd.iy(3),nlayer,length(hash(i0).index)];
    val=squeeze(...
        ncread(hash(i0).filename,field_name,ncstart,nccount)...
        );
  
    %Time index
    itstart=find(min(hash(i0).time)==t_roms,1,'first'); 
    it1=it-itstart+1; 
    it1(it1<1)=0; it1(it1>length(hash(i0).time))=0; it1(isnan(it1))=0; 
    
    
    %Interpolate
    for i1=1:length(out.val)
        for j2=1:size(wt,1)
            if it1(j2,i1)<1; continue; end
            for j1=1:size(w,1)
                w1=w(j1,i1)*wt(j2,i1); 
                if w1<1e-8; continue; end
                out.val(i1)=out.val(i1)+w1*...
                    val(ix(j1,i1),iy(j1,i1),it1(j2,i1));
                w_sum(i1)=w_sum(i1)+1; 
            end
        end
    end
    

end
out.val(w_sum<.99)=NaN; 

end

