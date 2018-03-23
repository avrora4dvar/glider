function [out] = get_roms_layer(roms_dir,field_name,loc,filter)
%GET_ROMS_LAYER Read layer from history files
%   [out]=get_roms_layer(roms_dir,field_name,loc,filter) reads layer of quantity field_name from
%   history files in the directory roms_dir. loc is a struct with the field
%   layer with the layer index and t with the output times [datenum].


if strcmpi(filter,'delta')
    tLim(1)=floor(min(loc.t)*24)/24;
    tLim(2)=ceil(max(loc.t)*24)/24;   
elseif strncmpi(filter,'mean',4)
    tPeriod=sscanf(filter,'mean%f')/24;
    tLim(1)=floor( (min(loc.t)-.5*tPeriod)*24 )/24;
    tLim(2)=ceil( (max(loc.t)+.5*tPeriod)*24 )/24;  
end
tLim=max(tLim,1/24); 
if isempty(tLim)
    error('time outside range'); 
end

%get roms indices
if false
    [hash,t_roms]=get_hash_roms(roms_dir,tLim);
else
    [hash,t_roms]=get_hash_roms(roms_dir,tLim,'ocean_avg');
end
    


%read roms
info=ncinfo(hash(1).filename,field_name); 
val1=nan(info.Size(1),info.Size(2),length(t_roms)); 
for i0=1:length(hash)
    iBnd(1)=find(t_roms==hash(i0).time(1),1); 
    iBnd(2)=find(t_roms==hash(i0).time(end),1); 
    if length(info.Size)==3
        val1(:,:,iBnd(1):iBnd(2))=squeeze(...
            ncread(hash(i0).filename,field_name,[1,1,hash(i0).index(1)],...
            [Inf,Inf,length(hash(i0).index)]) );
    elseif length(info.Size)==4
        val1(:,:,iBnd(1):iBnd(2))=squeeze(...
            ncread(hash(i0).filename,field_name,[1,1,loc.layer,hash(i0).index(1)],...
            [Inf,Inf,1,length(hash(i0).index)]) );
    end
end


%calculate weighted averages
if strcmpi(filter,'delta')
    val=nan(size(val1,1),size(val1,2),length(loc.t)); 
    wi=interp1(t_roms,[1:length(t_roms)],loc.t); 
   
    for it=1:length(loc.t)
        w=mod(wi(it),1); 
        val(:,:,it)=(1-w)*val1(:,:,floor(wi(it)))+...
            w*val1(:,:,ceil(wi(it))); 
    end
   
elseif strncmpi(filter,'mean',4)
    
    for it=1:length(loc.t)
        
        iBnd(1)=find(t_roms<=loc.t(it)-.5*tPeriod,1,'last');
        iBnd(2)=find(t_roms>=loc.t(it)+.5*tPeriod,1,'first');
        for i0=iBnd(1):iBnd(2)-1
            tmin=max(loc.t(it)-.5*tPeriod,t_roms(i0));
            tmax=min(loc.t(it)+.5*tPeriod,t_roms(i0+1));
            w(i0)=w(i0)+1/(t_roms(i0+1)-t_roms(i0))*...
                ( t_roms(i0+1)*(tmax-tmin)-.5*(tmax^2-tmin^2) );
            w(i0+1)=w(i0+1)+1/(t_roms(i0+1)-t_roms(i0))*...
                ( .5*(tmax^2-tmin^2)-(tmax-tmin)*t_roms(i0) );
        end
        w=w/sum(w);
        
        w=reshape(w(iBnd(1):iBnd(2)),[1,1,iBnd(2)-iBnd(1)+1]);
        w=repmat(w,[size(val1,1),size(val1,2),1]);
        val(:,:,it)=nansum(w.*val1(:,:,iBnd(1):iBnd(2)),3);
        
    end
   
end

out=struct('t',loc.t,'layer',loc.layer,'filter',filter,'field',field_name,'val',val); 

end

