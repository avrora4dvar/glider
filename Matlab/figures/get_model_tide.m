function [ zTide ] = get_model_tide(tide,loc)
%get_model_tide Calculates tide based on tidal consituents
%   [zTide]=get_model_tide(tide,loc); 

addpath('/home/server/student/homes/ipasmans/Matlab/roms/t_tide');

%initiate output
zTide=nan(size(loc.t)); 

if ~isempty(loc.t)
    
    %collect spatial points
    loc.lon=mod(loc.lon,360); loc.lon(loc.lon>180)=loc.lon(loc.lon>180)-360;
    xy=[loc.lon(:),loc.lat(:)];
    [xyUni,xy2uni,uni2xy]=unique(xy,'rows');
    
    for iCom=1:length(tide.freq)
        Afield=squeeze(tide.val(:,:,iCom)); Afield=permute(Afield,[2,1]);
        dfield=squeeze(tide.dval(:,:,iCom)); dfield=permute(dfield,[2,1]);
        A(:,iCom)=interp2(tide.lon',tide.lat',Afield,xyUni(:,1),xyUni(:,2));
        dA(:,iCom)=interp2(tide.lon',tide.lat',dfield,xyUni(:,1),xyUni(:,2));
    end
    
    for i1=1:size(xyUni,1)
        in=uni2xy==i1;
        tidecon=[abs(A(i1,:)); 0*abs(dA(i1,:));...
            rad2deg( angle(A(i1,:)) ); 0*rad2deg( angle(dA(i1,:)) )];
        tidecon=tidecon';
        
        if any(isnan(tidecon(:))); continue; end
        
        zTide(in)=t_predic(loc.t(in),tide.name,tide.freq,tidecon,'latitude',xyUni(i1,2));
    end
    
end

end

