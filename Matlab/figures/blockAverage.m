function [valOut] = blockAverage(varargin)
%BLOCKAVERAGE 
% valOut=blockAverage(tIn,valIn,tOut) averages valIn taking at tIn in
% blocks centered around tOut

%input
if length(varargin)==3
    tIn=varargin{1}; 
    valIn=varargin{2}; 
    tOut=varargin{3}; 
    if ~ismatrix(valIn)
        dim=find(size(valIn)~=1,1); 
    else
        dim=1; 
    end    
elseif length(varargin)==4
    tIn=varargin{1}; 
    valIn=varargin{2};
    tOut=varargin{3}; 
    dim=varargin{4}; 
end


%boundaries time blocks
tOut=reshape(tOut,[],1); 
tBnd=[-Inf; .5*tOut(1:end-1)+.5*tOut(2:end) ;Inf]; 

%permute dimensions
sIn=size(valIn); 
dimIndex=[1:ndims(valIn)]; 
dimIndex=circshift(dimIndex,[0,1-dim]); 
valIn=permute(valIn,dimIndex); 
sTmp=size(valIn);
valIn=reshape(valIn,size(valIn,1),[]);
maskIn=ones(size(valIn)); maskIn(isnan(valIn))=0; 

%initiate output
sOut=size(valIn); sOut(1)=length(tOut); 
valOut=zeros(sOut); 

%calculate block average
for i1=1:length(tBnd)-1
   maskOut=zeros(1,size(valOut,2));
    
   inBlock= tIn==tBnd(i1) | tIn==tBnd(i1+1);
   valOut(i1,:)=valOut(i1,:)+.5*nansum(valIn(inBlock,:),1); 
   maskOut=maskOut+.5*sum(maskIn(inBlock,:),1); 
   
   inBlock= tIn>tBnd(i1) & tIn<tBnd(i1+1); 
   valOut(i1,:)=valOut(i1,:)+nansum(valIn(inBlock,:),1); 
   maskOut=maskOut+sum(maskIn(inBlock,:),1); 

   valOut(i1,:)=valOut(i1,:)./maskOut;
   valOut(i1,maskOut==0)=NaN; 
end

%convert output back
sTmp(1)=size(valOut,1); 
valOut=reshape(valOut,sTmp); 
valOut=ipermute(valOut,dimIndex); 


end

