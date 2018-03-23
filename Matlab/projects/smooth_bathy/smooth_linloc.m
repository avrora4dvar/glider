function [h] = smooth_lin(h0,mask,rlim)
%smooth_meo smooths bathymetry using Mellor-Ezer-Oey scheme (Mellor 1994). 
%  h=smooth_meo(h0,mask,r) gives a smooth bathymetry h from a non-smooth
%  bathymetry h0, with mask the mask of active points and r the smoothing
%  parameter (r=0->uniform bathymetry). 

nlayer=40; 

warning off all
opt=optimset('display','off'); 

s=size(h0); 
h=h0; 
nin=sum(mask(:)); 

[ix,iy]=ndgrid([1:size(h,1)],[1:size(h,2)]); ix=ix(:); iy=iy(:); 
j1list=[0,0,1,-1]; j2list=[1,-1,0,0]; 

iterNo=0;  r=1; 
while true
    iterNo=iterNo+1; 
    
    %create matrices
    Id=eye(9);
    A=zeros(12*2,9);
    ind=[1:9]; ind=reshape(ind,[3,3]);
    j0=1;
    for j1=1:3
        for j2=1:2
            A(j0,ind(j1,j2))=-r-1;
            A(j0,ind(j1,j2+1))=-r+1;
            A(j0+1,ind(j1,j2))=-r+1;
            A(j0+1,ind(j1,j2+1))=-r-1;
            j0=j0+2;
        end
    end
    for j2=1:3
        for j1=1:2
            A(j0,ind(j1,j2))=-r-1;
            A(j0,ind(j1+1,j2))=-r+1;
            A(j0+1,ind(j1,j2))=-r+1;
            A(j0+1,ind(j1+1,j2))=-r-1;
            j0=j0+2;
        end
    end
    
    
    %Beckman_Haidvogel number
    rx0=beckman_haidvogel_haney(h);
    maxr=max(rx0(:));
    [i1,i2]=find(rx0>r); 
    
    %reshuffle
    w=randperm(size(i1,1)); 
    i1=i1(w); i2=i2(w);     
    
    for i0=1:size(i1,1)
        
        hloc=nan(9,1); 
        hfix=zeros(9,1);
        for k1=1:3
            for k2=1:3
                j1=i1(i0)+k1-2; j2=i2(i0)+k2-2; 
                if j1<1||j2<1||j1>s(1)||j2>s(2)
                    hloc(ind(k1,k2))=NaN;
                elseif ~mask(j1,j2)
                    hloc(ind(k1,k2))=NaN; 
%                     hloc(ind(k1,k2))=h(j1,j2); 
%                     hfix(ind(k1,k2))=1; 
                else
                    hloc(ind(k1,k2))=h(j1,j2);
                end
                
            end
        end
        hnan=~isnan(hloc); 
        hfix=hloc<50; 
         
        inan=zeros(size(A,1),1); inan=false; 
        for j1=1:size(A,1)
            inan(j1)=~any( A(j1,:)~=0 & isnan(hloc') ); 
        end
        hloc(isnan(hloc))=0; 
        if ~any(inan); continue; end
        

        %current slopes
        Aloc=A(inan,:); 
        Iloc=exp( (50./hloc).^2 ); Iloc=max(Iloc,1); Iloc(~hnan)=1;
        Iloc=diag(Iloc); 
        Zloc=zeros(size(Aloc,2),1); 
        bloc=-Aloc*hloc;   
        
        %add fixes
%         for k0=1:length(hfix)
%            if hfix(k0)
%               Aloc=[Aloc;zeros(2,size(Aloc,2))]; 
%               Aloc(end-1,k0)=1; 
%               Aloc(end,k0)=-1;
%               bloc=[bloc;0;0];
%            end
%         end
%         
        
        %minimize
        dh=quadprog(Iloc,Zloc,Aloc,bloc,[],[],[],[],Zloc,opt); 
        hout=hloc+dh; 
        
        for k1=1:3
            for k2=1:3
                j1=i1(i0)+k1-2; j2=i2(i0)+k2-2;
                if j1<1||j2<1||j1>s(1)||j2>s(2)
                    continue; 
                elseif ~mask(j1,j2)
                    continue; 
                else
                    h(j1,j2)=hout(ind(k1,k2));
                end
            end
        end

    end
    
    if mod(iterNo,1)==0
        %Calculate Haney number
        sigma=[-nlayer+.5:-.5]/nlayer;
        z=roms_sigma2z(sigma,h); 
        rx0=beckman_haidvogel_haney(h); maxr=max(rx0(:)); 
        rx1=beckman_haidvogel_haney(z); maxr1=prctile(rx1(:),99.95); 
        
        display(sprintf('Iteration %d. Max. rx0=%e Max. rx1=%e',iterNo,maxr,maxr1 ));
        if maxr1<rlim; break; end
    end

    %update r
    r=r*.975; 
    

end


end

