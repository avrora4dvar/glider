function [ out ] = get_roms_3d(varargin)
%GET_ROMS_3D
% out=get_roms_3d(grid_file,roms_dir,field_name,grd,loc) reads field
% field_name from history files in roms_dir at locations in loc
% Here grd is struct with fields lon,lat,mask containing the grid
% associated with field_name and loc a struct with fields lon,lat,t,z/depth
% containing the longitudes,latitudes,times (in datenum) and z/depth below
% surface where field_name is requested. 

%input
grid_file=varargin{1}; 
roms_dir=varargin{2}; 
field_name=varargin{3}; 
grd=varargin{4}; 
loc=varargin{5}; 

npoints=length(loc.lon);

%get grid indices
[ix,iy,iz,w,loc]=loc2index_3d(grid_file,roms_dir,grd,loc); 
iBnd.ix(1)=min(ix(:)); iBnd.ix(2)=max(ix(:)); 
iBnd.ix(3)=iBnd.ix(2)-iBnd.ix(1)+1; 
iBnd.iy(1)=min(iy(:)); iBnd.iy(2)=max(iy(:)); 
iBnd.iy(3)=iBnd.iy(2)-iBnd.iy(1)+1;
iBnd.iz(1)=min(iz(:)); iBnd.iz(2)=max(iz(:)); 
iBnd.iz(3)=iBnd.iz(2)-iBnd.iz(1)+1; 

%Shift domain
ix=ix-iBnd.ix(1)+1;
iy=iy-iBnd.iy(1)+1; 
iz=iz-iBnd.iz(1)+1; 
ix(isnan(ix))=1; iy(isnan(iy))=1; iz(isnan(iz))=1; 

%get roms indices
[hash,t_roms]=get_hash_roms(roms_dir,[min(loc.t),max(loc.t)],'ocean_his');

%Time index
[it,wt]=loc2index_1d(t_roms,loc.t);

%read data
out=struct('lon',loc.lon,'lat',loc.lat,'z',loc.z,'depth',loc.depth,'t',loc.t,'field',field_name,'val',zeros(size(loc.t))); 
w_sum=zeros(size(loc.t));
for i0=1:length(hash)
    ncstart=[iBnd.ix(1),iBnd.iy(1),iBnd.iz(1),hash(i0).index(1)];
    nccount=[iBnd.ix(3),iBnd.iy(3),iBnd.iz(3),length(hash(i0).index)];
    val=squeeze(...
        ncread(hash(i0).filename,field_name,ncstart,nccount)...
        );
    
    %Time index
    itstart=find(hash(i0).time(1)==t_roms,1,'first');
    it1=it-itstart+1;
    it1(it1<1)=0; it1(it1>length(hash(i0).time))=0; it1(isnan(it1))=0;
    
    %Interpolate
    for i1=1:length(out.val)
        if isnan(it(1,i1)); continue; end
        for j2=1:size(wt,1)
            if it1(j2,i1)<1; continue; end
            for j1=1:size(w,1)
                w1=w(j1,i1)*wt(j2,i1); 
                if w1<1e-8; continue; end
                out.val(i1)=out.val(i1)+w1*...
                    val(ix(j1,i1),iy(j1,i1),iz(j1,i1),it1(j2,i1));
                w_sum(i1)=w_sum(i1)+w1; 
            end
        end
    end 
    
end
out.val(w_sum<.99)=NaN;        

end

