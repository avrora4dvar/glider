function [valOut,cIn] = filter_running_tide(tIn,valIn,dim)
%FILTER_RUNNING_MEAN valOut=filter_running_mean(tIn,valIn,dim)
%   Running mean filter for regular grid averaging over [-ndt:ndt] points

%Tides M2,S2,K1,O1
omega=deg2rad([28.984104,30,15.041069,13.943035]*24); 
%Tides M2-S2, K1-O1
omega=[mean(omega(1:2)),mean(omega(3:4))]; 

%reshape input
[valIn,valSize]=array2matrix(valIn,dim); 
tIn=reshape(tIn,[],1); 

%Create tidal matrix
M=ones(length(tIn),2+2*length(omega)); 
M(:,2)=tIn(:)-.5*(min(tIn)+max(tIn));  
if all(imag(valIn(:))==0)
    display('Scalar tidal fit'); 
    for i2=1:length(omega)
        M(:,2*i2+1)=cos(omega(i2)*tIn); M(:,2*i2+2)=sin(omega(i2)*tIn); 
    end
else
    display('Field tidal fit'); 
    for i2=1:length(omega)
        M(:,2*i2+1)=exp(1i*omega(i2)*tIn); M(:,2*i2+2)=exp(-1i*omega(i2)*tIn); 
    end
end

%Times for which tidal fit will be calculated
dt=max(360./omega)/24; 
tOut=[min(tIn)+.5*dt:dt:max(tIn)-.5*dt]; 


%Fit tides
c=nan(size(M,2),length(tOut),size(valIn,2)); 
for i1=1:size(c,2)
    for i2=1:size(c,3)
        in=~isnan(valIn(:,i2)) & tIn>=tOut(i1)-dt & tIn<=tOut(i1)+dt; 
        if sum(in)<4*size(M,2); continue; end
        c(1:end,i1,i2)=M(in,1:end)\valIn(in,i2);  
    end
end

%Interp c
c=permute(c,[2,1,3]);
cIn=nan(length(tIn),size(c,2),size(c,3)); 
for i3=1:size(c,3)
    in=~isnan(c(:,1,i3));
    if ~any(in); continue; end
    cIn(:,:,i3)=interp1(tOut(in),c(in,:,i3),tIn); 
    iL=find(in,1,'first'); iU=find(in,1,'last');  
    for i2=1:size(c,2)
     cIn(tIn<=tOut(iL),i2,i3)=c(iL,i2,i3); 
     cIn(tIn>=tOut(iU),i2,i3)=c(iU,i2,i3);
    end
end
cIn=permute(cIn,[2,3,1]); 

%Generate tide signal
valOut=nan(size(valIn)); 
for i3=1:size(cIn,3)
    valOut(i3,:)=M(i3,3:end)*cIn(3:end,:,i3); 
end

%Remove tide
valOut=valIn-valOut; 

%Reshape back
valOut=matrix2array(valOut,dim,valSize); 
    
end

