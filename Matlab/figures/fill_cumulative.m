function [valOut,valCum] = fill_cumulative(tIn,valIn,tOut,dim,varargin)
%FILL_CUMULATIVE [valOut,valCum]=fill_cumulative(tIn,valIn,tOut)
%   Create piecewise linear interpolation of with the same piecewise
%   average as the original signal

iPar=find(strcmpi('average',varargin));
if ~isempty(iPar)
    flag_avg=true; 
else
    flag_avg=false; 
end

%Make dim the 1st dimension
[valIn,valSize]=array2matrix(valIn,dim); 
tOut=reshape(tOut,[],1); 

%Initiate output 
valOut=nan(length(tOut),size(valIn,2));
valCum=nan(length(tOut)+2,size(valIn,2)); 
valSize(dim)=length(tOut); 

%Numerical integration
dt=diff(tOut); dt=[dt;0]; 
for i2=1:size(valOut,2) 
   
   if ~flag_avg
       in=~isnan(valIn(:,i2));
       Qcum=cumtrapz(tIn(in),valIn(in,i2));
   else
       in=find(~isnan(valIn(:,i2))); 
       Qcum=.5*diff(tIn(in)).*valIn(in(1:end-1),i2)...
           +.5*diff(tIn(in)).*valIn(in(2:end),i2); 
       Qcum=[0;cumsum(Qcum)]; %add start time
   end
   
   itL=find(tOut>=min(tIn(in)),1,'first'); itU=find(tOut<=max(tIn(in)),1,'last'); 
   if isempty(itL) || isempty(itU); continue; end

   valCum(itL+1:itU+1,i2)=interp1(tIn(in),Qcum,tOut(itL:itU));
   valCum([itL,itU+2],i2)=interp1(tIn(in),valIn(in,i2),[tOut(itL),tOut(itU)]); 
   
   %Create spline
   pp=spline(tOut(itL:itU),valCum(itL:itU+2,i2));
   pp.coefs=bsxfun(@times,pp.coefs,[3,2,1,0]); 
   pp.coefs=circshift(pp.coefs,[0,1]); 
   
   %Get values from spline
   valOut(itL:itU,i2)=ppval(pp,tOut(itL:itU));  
   valOut(1:itL,i2)=valOut(itL,i2); 
   valOut(itU:end,i2)=valOut(itU,i2); 
end

%Transform dimension back
valCum=valCum(2:end-1,:);
valOut=matrix2array(valOut,dim,valSize); 
valCum=matrix2array(valCum,dim,valSize); 
    
end

