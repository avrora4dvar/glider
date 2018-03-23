function [ h ] = smooth_ext(h,rmax)
%smooth_max Summary of this function goes here
%   Detailed explanation goes here

dt=.01; dx=.005;  nanr=.5;

rx=beckman_haidvogel_haney(h); 
h0=h; 

niter=0; 
while max(rx(:))>rmax
    display(sprintf('rmax:%.2f',max(rx(:)))); 
    
    if mod(niter,20)==0 && niter>0
       [ix,iy]=find(rx>rmax); 
       for k=1:length(ix)
          if ix(k)==1 || ix(k)==size(h,1) || iy(k)==1 || iy(k)==size(h,2); continue; end
          val1=[h(ix(k)-1,iy(k)),h(ix(k)+1,iy(k)),h(ix(k),iy(k)-1),h(ix(k),iy(k)+1)]; 
          if any(isnan(val1)); h(ix(k),iy(k))=NaN; end
       end
       [ix,iy]=find(~isnan(h)); 
       for k=1:length(ix)
          if ix(k)==1 || ix(k)==size(h,1) || iy(k)==1 || iy(k)==size(h,2); continue; end
          val1=[h(ix(k)-1,iy(k)),h(ix(k)+1,iy(k)),h(ix(k),iy(k)-1),h(ix(k),iy(k)+1)]; 
          if sum(isnan(val1))>=3; h(ix(k),iy(k))=NaN; end
       end
       [ix,iy]=find(isnan(h)); 
       for k=1:length(ix)
          if ix(k)==1 || ix(k)==size(h,1) || iy(k)==1 || iy(k)==size(h,2); continue; end
          val1=[h(ix(k)-1,iy(k)),h(ix(k)+1,iy(k)),h(ix(k),iy(k)-1),h(ix(k),iy(k)+1)]; 
          if sum(~isnan(val1))>=3; h(ix(k),iy(k))=nanmean(val1); end
       end
    end
    
    %lon flux
    mhu=.5*h(1:end-1,:)+.5*h(2:end,:);
    dhu=h(2:end,:)-h(1:end-1,:);
    fluxu=.5*( rmax*mhu-dhu )*dt;
    fluxu( abs(dhu./mhu)<rmax )=0;
    fluxu=sign(fluxu).*min(abs(fluxu),nanmin(h(1:end-1,:),h(2:end,:))*dx); 
    fluxu(isnan(fluxu))=0; 
   

    %lat flux
    mhv=.5*h(:,2:end)+.5*h(:,1:end-1);
    dhv=h(:,2:end)-h(:,1:end-1);
    fluxv=.5*( rmax*mhv-dhv )*dt;
    fluxv( abs(dhv./mhv)<rmax  )=0;
    fluxv=sign(fluxv).*min(abs(fluxv),nanmin(h(:,1:end-1),h(:,2:end))*dx);
    fluxv(isnan(fluxv))=0; 
    
    
    %update
    h(2:end,:)=h(2:end,:)+fluxu;
    h(1:end-1,:)=h(1:end-1,:)-fluxu;
    h(:,2:end)=h(:,2:end)+fluxv;
    h(:,1:end-1)=h(:,1:end-1)-fluxv;
    
    
    rx=beckman_haidvogel_haney(h); 
    niter=niter+1; 
        
    
end

display(sprintf('N iter: %d',niter)); 


    
end

