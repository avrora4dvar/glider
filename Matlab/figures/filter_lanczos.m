function yf = filter_lanczos(y,dt,fcut,order,dim)
%FILTER_LANCZOS LP-filter data using cos-Lanczos filter
%[yf]=filter_lanczos(y,dt,fcut,order,dim) filters data in y along dimension
%dim using a window that is 2*order+1 points wide. fcut is the
%half-amplitude frequency and dt the time step in y. 

%Stencil
t=[-order:order]'; 
fcut=2*fcut; 
w=0.5*(1+cos(pi*t/order)).*sinc(dt*t*fcut);
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


% %convolution signal with weights
% yf=nan(size(y));
% y=permute(y,[2,1]); yf=permute(yf,[2,1]); ynan=permute(ynan,[2,1]); 
% for i2f=1:size(yf,2)
%    iBndY(1)=max(1,i2f-order); 
%    iBndY(2)=min(size(y,2),i2f+order); 
%    iBndW(1)=order+1-(i2f-iBndY(1)); 
%    iBndW(2)=order+1+(iBndY(2)-i2f);
%   
%    
%    if iBndW(2)>=iBndW(1)
%        yf(:,i2f)=0; 
%        wf=zeros(size(yf(:,i2f))); 
%        for i2=0:iBndW(2)-iBndW(1)
%            yf(:,i2f)=yf(:,i2f)+...
%                w(i2+iBndW(1))*y(:,i2+iBndY(1)); 
%            wf=wf+...
%                w(i2+iBndW(1))*ynan(:,i2+iBndY(1)); 
%        end
%        yf(:,i2f)=yf(:,i2f)./wf;
%        
%        yf(sum(ynan(:,iBndY(1):iBndY(2)),2)<order*1.7,i2f)=NaN; 
%    end
% end
% yf=ipermute(yf,[2,1]); clear y; clear ynan;  

%reshape back
yf=matrix2array(yf,dim,s); 

end

