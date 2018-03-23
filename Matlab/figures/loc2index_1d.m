function [it,w] = loc2index_1d(tBase,tIn)
%LOC2INDEX [it,w]=loc2index_1d(tBase,tIn)
%   Indices and weights 1d array

%Index 
it1=interp1(reshape(tBase,1,[]),[1:length(tBase)],reshape(tIn,1,[])); 
it(1,:)=floor(it1); 
it(2,:)=ceil(it1); 
    
%Weight
w(1,:)=1-abs(it(1,:)-it1); w(2,:)=1-abs(it(2,:)-it1); 
w(isnan(w))=0; 

%Normalize
wsum=sum(w,1);
in=wsum>0; 
w(:,in)=bsxfun(@times,w(:,in),1./wsum(in)); 
    
    
    
end

