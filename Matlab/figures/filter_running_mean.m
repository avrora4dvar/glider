function yf = filter_running_mean(y,order,dim)
addpath('/home/server/student/homes/ipasmans/Matlab/figures'); 

%FILTER_LANCZOS LP-filter data using cos-Lanczos filter
%[yf]=filter_lanczos(y,dt,fcut,order,dim) filters data in y along dimension
%dim using a window that is 2*order+1 points wide. fcut is the
%half-amplitude frequency and dt the time step in y. 

%Stencil
t=[-order:order]'; 
w=ones(size(t)); 
w=reshape(w/sum(w),[],1); %normalize

%reshape
[y,s]=array2matrix(y,dim); 

%nan
ynan=ones(size(y)); 
ynan(isnan(y))=0; 
y(isnan(y))=0; 

%Output
yf=conv2(y,w,'same'); 
ynan=conv2(ynan,w,'same'); 
yf=yf./ynan; 
yf(ynan<.5)=NaN; 

%reshape back
yf=matrix2array(yf,dim,s); 

end

%% Array <-> matrix

function [y,s]=array2matrix(y,dim) 
    diml=[1:ndims(y)]; 
    diml=circshift(diml,[0,-dim+1]); 
    y=permute(y,diml); 
    s=size(y); 
    y=reshape(y,size(y,1),[]); 
end

function y=matrix2array(y,dim,s)
   s(1)=size(y,1); 
   y=reshape(y,s); 
   diml=[1:ndims(y)]; 
   diml=circshift(diml,[0,-dim+1]); 
   y=ipermute(y,diml);    
end
    

