function [h] = smooth_mb(h0,mask,r,alpha)
%smooth_mb smooths bathymetry using Matinho&Batten (2006) scheme. 
%  h=smooth_mb(h0,mask,r) gives a smooth bathymetry h from a non-smooth
%  bathymetry h0, with mask the mask of active points and r the smoothing
%  parameter (r=0->bathymetry is left unaltered). 

s=size(h); 
h=h0; 
h(~mask)=NaN; 
nin=sum(mask(:)); 

[ix,iy]=ndgrid([1:size(h,1)],[1:size(h,2)]); ix=ix(:); iy=iy(:); 

%Beckman_Haidvogel number
rx0=beckman_haidvogel_haney(h); 

iterNo=1; 
while (sum(rx0(:)>r)/nin>alpha)
    
    maxr=max(rx0(:)); 
    display(sprintf('iteration %d. Max. rx',iterNo,maxr); 
    iterNo=iterNo+1; 
    
    [i1,i2]=find(rx0==maxr,1); 
    
    
    for j1=[i1-1,i1+1]
        if j1<1||j1>s(1); continue; end
        if ~ismask(j1,i2); continue; end
        rloc=-(h(j1,i2)-h(i1,i2))/(h(i1,i2)+h(j1,i2));
        if rloc==rmax
            h(j1,i2)=h(i1,i2)*(1-r)/(1+r); 
        elseif rloc==-rmax
            h(i1,i2)=h(j1,i2)*(1+r)/(1-r); 
        end    
    end
    
    for j2=[i2-1,i2+1]
        if j2<1||j2>s(2); continue; end
        if ~ismask(i1,j2); continue; end
        rloc=-(h(i1,j2)-h(i1,i2))/(h(i1,i2)+h(i1,j2));
        if rloc==rmax
            h(i1,j2)=h(i1,i2)*(1-r)/(1+r);
        elseif rloc==-rmax
            h(i1,i2)=h(i1,j2)*(1+r)/(1-r);
        end
    end
        
    %Beckman_Haidvogel number
    rx0=beckman_haidvogel_haney(h);
end


end

