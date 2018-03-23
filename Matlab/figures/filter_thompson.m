function [y,R,w]=filter_thompson(x,nDim,filterData)
%DETIDER low-pass tidal filter based on Thompson (1983)
%
%Syntax: 
% y=detider(x,n,filterData)
%Input:
% x: [nx1 DBL] input time-series
% n: [INT] order of filter. n-order filter uses 2xn+1 filter coefficients
% filterData: struct with fields
% -f_s: sample frequency
% -f_pass: frequency below which filter respons approximates 1
% -f_stop: frequency above which filter respons approximates 0
% -f [1xm DBL,optional]: tidal frequencies
% -A [1xm DBL,optional]: tidal amplitudes. Do only use if f is specified
% -df [DBL/1xm DBL,optional]: half-width frequency range associated with
%  the tidal frequencies
%Frequencies must be specified in the same units as f_s
%Output:
% y=filtered data
% R=struct with fields
%  -f: [1xl DBL] normalized frequency (with 1 0.5*f_s)
%  -H: [1xl DBL] filter transfer for given frequency
% w=[2xn+1 DBL] filter coefficients
%Example:
% Example of hourly samplied time series
% filterData=struct('fs',1/24,'f_pass',1/36,'f_stop',1/25,'A',1,'f',1/12,'df',.01);
% t=[0:100]/24; 
% x=1*sin(2*pi*t/12)+.1*cos(2*pi*t/4); 
% y=detider(x,24,filterData); 

%check input
if filterData.f_pass>filterData.f_stop
    error('f_pass must be smaller than f_stop'); 
end
if length(filterData.f)~=length(filterData.A)
    error('number of tidal frequencies and amplitudes must be equal'); 
end
n=filterData.order;
if isempty(n) || n<1
    error('filter order must be larger than zero'); 
end

%normalize frequencies
omega=2*pi*filterData.f(:)/filterData.f_s;
dOmega=2*pi*filterData.df(:);
if length(dOmega)==1
    dOmega=dOmega*ones(size(omega)); 
end
omega2=2*pi*filterData.f_stop/filterData.f_s;
omega1=2*pi*filterData.f_pass/filterData.f_s; 

%weighting coefficients
if isempty(omega)
    a=1;
    b=1;
else
    a= 1*sqrt(sum(filterData.A.^2));
    b= 1*sqrt(sum(filterData.A.^2));
end

%Hessian
H=a^2*(intR2mat(omega1,n)-intR2mat(0,n))/omega1+...
    b^2*(intR2mat(pi,n)-intR2mat(omega2,n))/(pi-omega2);
for i1=1:length(omega)
    H=H+filterData.A(i1)^2*(intR2mat(omega(i1)+dOmega(i1),n)-intR2mat(omega(i1)-dOmega(i1),n))/(2*dOmega(i1)); 
end

%gradient
g=2*a^2*(intRmat(omega1,n)-intRmat(0,n)); 
g=g'; 

%weights
H=real(H); g=real(g); 
w=H\g;  
w(1)=2*w(1); 
w=[flipdim(w(2:end),1);w]; 
w=w/sum(w); 

%transfer function  
R.f=[0:.001:.099,.1:.01:1]; 
for i1=1:length(R.f)
        R.H(i1)=sum( w'.*exp(2i*pi*R.f(i1)*.5*[-n:n]) ); 
end
R.H=real(R.H); 

%filtered output
y=nan(size(x)); 
dimX=[1:length(size(x))];  dimPer=circshift(dimX,[0,-nDim+1]); 
x=permute(x,dimPer); y=permute(y,dimPer);  
dimShape=size(x); 
x=reshape(x,size(x,1),[]); y=reshape(y,size(y,1),[]); 
w=repmat(w,[1,size(x,2)]); 
for i1=n+1:size(x,1)-n
    y(i1,:)=sum(w.*x(i1-n:i1+n,:),1); 
end
% for i1=1:n
%     ww=w(nDim+1-(i1-1):end); ww=ww/sum(ww); 
%     y(i1,:)=sum(x(1:length(ww)).*ww);
% end
x=reshape(x,dimShape); y=reshape(y,dimShape); 
x=ipermute(x,dimPer); y=ipermute(y,dimPer); 


end

%% support functions

function y=Rmat(w,n)
%Rmat*w=transfer function

y=exp(1i*[0:n]*w)+exp(-1i*[0:n]*w); 
end

function y=intRmat(w,n)
%Frequency primitive of Rmat

y=[2*w,exp(1i*[1:n]*w)./(1i*[1:n])+exp(-1i*[1:n]*w)./(-1i*[1:n])]; 
end

function y=intR2mat(w,n)
%Frequency primitive of Rmat'*Rmat

    [expMat1,expMat2]=ndgrid([0:n],[0:n]);
    
    expMat=nan(n+1,n+1,4); 
    expMat(:,:,1)=-expMat1+expMat2;
    expMat(:,:,2)=-expMat1-expMat2; 
    expMat(:,:,3)=expMat1+expMat2; 
    expMat(:,:,4)=expMat1-expMat2; 
    expMat=1i*expMat; 
    
    y=exp(w*expMat)./expMat; 
    y(expMat==0)=w; 
    y=sum(y,3); 
end

