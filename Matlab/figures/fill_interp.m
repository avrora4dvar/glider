function [valOut] = fill_interp(tIn,valIn,tOut,dtMax,dim)
%FILL_INTERP [valOut]=fill_interp(tIn,valIn,tOut,dtMax,dim)
%  
addpath('/home/server/student/homes/ipasmans/Matlab/figures/'); 


%Make dim the 1st dimension
[valIn,valSize]=array2matrix(valIn,dim); 

%Initiate output
valOut=nan(length(tOut),size(valIn,2)); 
valSize(dim)=length(tOut); 


%Interpolate
for i2=1:size(valIn,2)
    %Interpolate
    in=~isnan(valIn(:,i2));
    if sum(in)<2
        valOut(:,i2)=NaN;
    else
        valOut(:,i2)=interp1(tIn(in),valIn(in,i2),tOut,'linear',NaN);
        dt=abs(tOut-interp1(tIn(in),tIn(in),tOut,'nearest',NaN));
        valOut(dt>dtMax | isnan(dt),i2)=NaN;
    end
end

%Back transformation
valOut=matrix2array(valOut,dim,valSize); 

end

