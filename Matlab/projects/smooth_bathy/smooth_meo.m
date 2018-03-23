function [h] = smooth_meo(h0,mask,rlim)
%smooth_meo smooths bathymetry using Mellor-Ezer-Oey scheme (Mellor 1994). 
%  h=smooth_meo(h0,mask,r) gives a smooth bathymetry h from a non-smooth
%  bathymetry h0, with mask the mask of active points and r the smoothing
%  parameter (r=0->uniform bathymetry). 

nlayer=40; 

s=size(h0); 
h=h0; 
h(~mask)=NaN; 
nin=sum(mask(:)); 

[ix,iy]=ndgrid([1:size(h,1)],[1:size(h,2)]); ix=ix(:); iy=iy(:); 
j1list=[0,0,1,-1]; j2list=[1,-1,0,0]; 


iterNo=0;  r=1; 
while true
    iterNo=iterNo+1; 
    
    %Beckman_Haidvogel number
    rx0=beckman_haidvogel_haney(h);
    maxr=max(rx0(:));
    [i1,i2]=find(rx0>r); 
    
    %reshuffle
    w=randperm(size(i1,1)); 
    i1=i1(w); i2=i2(w);     
    w=interp1([0,1],[.5,4+.5],rand(size(i1,1),1));  w=round(w);
    
    for i0=1:size(i1,1)
        j1list=circshift(j1list,[0,w(i0)]); 
        j2list=circshift(j2list,[0,w(i0)]); 
        for j0=1:4
            
            %skip conditions
            j1=i1(i0)+j1list(j0); j2=i2(i0)+j2list(j0);
            if j1<1 || j2<1 || j1>s(1) || j2>s(2); continue; end
            if ~mask(i1(i0),i2(i0))||~mask(j1,j2); continue; end
            
            %smooth
            rloc=(h(i1(i0),i2(i0))-h(j1,j2))/(h(i1(i0),i2(i0))+h(j1,j2));
            if rloc>r
                delta=.5*(h(i1(i0),i2(i0))-h(j1,j2))-.5*r*(h(i1(i0),i2(i0))+h(j1,j2));
                h(i1(i0),i2(i0))=h(i1(i0),i2(i0))-delta;
                h(j1,j2)=h(j1,j2)+delta;
            elseif rloc<-r
                delta=.5*(h(j1,j2)-h(i1(i0),i2(i0)))-.5*r*(h(i1(i0),i2(i0))+h(j1,j2));
                h(j1,j2)=h(j1,j2)-delta;
                h(i1(i0),i2(i0))=h(i1(i0),i2(i0))+delta;
            end
            
        end
    end
    
    if mod(iterNo,10)==0
        %Calculate Haney number
        sigma=[-nlayer+.5:-.5]/nlayer;
        z=roms_sigma2z(sigma,h); 
        rx0=beckman_haidvogel_haney(h); maxr=max(rx0(:)); 
        rx1=beckman_haidvogel_haney(z); maxr1=prctile(rx1(:),99.98); 
        display(sprintf('Iteration %d. Max. rx0=%e Max. rx1=%e',iterNo,maxr,maxr1 ));
        if maxr<rlim; break; end
    end

    %update r
    r=r*.975; 
    

end


end

