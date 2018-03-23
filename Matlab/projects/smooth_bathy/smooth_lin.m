function x=smooth_lin(h,x,mask)


%% create input speesest gradient algorithm
W=ones(numel(h),1)/sum(mask(:)); 
%W=W.*max(1,(50./h(:)).^2); 

%settings for steepest gradient
b=h(:); 
ind=[1:length(b)]; ind=reshape(ind,size(h,1),size(h,2)); 

x=x(:); 
%x=mean(b)*ones(size(b)); 
%x=b; 

%% steepest gradient algorithm

iIn0=geoMask(ind,mask); 

%minimize J=(b-x)'*W*(b-x)-1/t sum_i log(-(Gx)_i)

Gx=G(ind,'N',x,iIn0); Gx=Gx(iIn0);
while any(Gx>=0)
	x=min(b)*ones(size(b))+0*(rand(size(x))-.5); 
    Gx=G(ind,'N',x,iIn0); Gx=Gx(iIn0); 
end
t=100000; 
eps=1e-100;

k=0;  R=[]; J1=[]; J2=[]; 
figure(); 
while true
    r=b-x; %residual
    
    R=[R,norm(r)]; 
    J1=[J1,sum(r.*W.*r)];
    J2=[J2,sum(log(max(-Gx,eps)))/t]; 
    
    %direction descent
    lam=G(ind,'N',x,iIn0); lam=-1./lam;  %1./Gx
    dir1=2*W.*r; dir2=G(ind,'T',lam,iIn0)/t; 
    dir=dir1-dir2; %-DJ(x)  
    
    %Constraint
    dWr=sum(dir.*W.*r); dWd=sum(dir.*W.*dir); 
    Gd=G(ind,'N',dir,iIn0);  Gd=Gd(iIn0); 
    Gx=G(ind,'N',x,iIn0);  Gx=Gx(iIn0);   

    
    %maximum gamma that keeps point in interior
    gamma=2*dWr/dWd;
    for i1=1:length(Gx)
        if -Gx(i1)-gamma*Gd(i1)<=0
            gamma=min(gamma,-Gx(i1)./Gd(i1)); 
        end
    end

   
    %find gamma for which J(x+gamma*dir) is minimal
    gamma=[0:.01:.99]*gamma;
    f=zeros(size(gamma));
    for i1=1:length(f)
        f(i1)=-2*gamma(i1)*dWr+gamma(i1)^2*dWd-sum(log( max(eps, -Gx-gamma(i1)*Gd) ),1)/t;
    end
    [dummy,iMin]=min(f); gamma=gamma(iMin);
    
    if mod(k,20)==0
        display(sprintf('Iteration %d. J1: %e. J2: %e.',k,J1(end),J2(end)  ));
        display(sprintf('Iteration %d. Norm: %e.',k,norm(r) ));
        display(sprintf('Iteration %d. Direction: %e %e',k,norm(dir1),norm(dir)));
        display(sprintf('Iteration %d. Gamma:%e %d',k,gamma,iMin));
        
        %Beckman_Haidvogel
        x1=reshape(x,size(h)); x1(~mask)=NaN;
        rx0=beckman_haidvogel_haney(x1);
        rx0=rx0(mask(:));
        display(sprintf('Iteration %d. Max rx: %e. Mean rx: %e',k,...
            max(rx0),mean(rx0) ));
    end
    
    %next x
    x=x+gamma*dir;  
    
    %iteration number
    k=k+1; 
    if k==200000; break; end
    if gamma<eps; break; end
 
end

x=reshape(x,size(h)); 


end


%%

function [iIn]=geoMask(ind,mask)
%[iIn]=geoMask(ind) find rows of G that connect 2 cells that lie both
%within the mask

j1list=[-1,1,0,0];
j2list=[0,0,-1,1];

iIn=ones(4,numel(ind));
for i2=1:size(ind,2)
    for i1=1:size(ind,1)
        
        for j0=1:length(j1list)
            j1=i1+j1list(j0); j2=i2+j2list(j0); 
            if j1<1 || j2<1 || j1>size(ind,1) || j2>size(ind,2)
                iIn(j0,ind(i1,i2))=0;
            elseif ~mask(i1,i2) || ~mask(j1,j2)
                iIn(j0,ind(i1,i2))=0;                 
            end
        end
        
    end
end

iIn=iIn(:)==1; 

end

%%

function [Gx]=G(ind,side,x,iIn)
%Gx=G(ind,side,x,iIn) Calculates G*x (side='N') or G'*x (side='T')

j1list=[-1,1,0,0]; 
j2list=[0,0,-1,1]; 
r=.18; %maximum Beckman-Haidvogel number
r=[-r-1,-r+1]; r=r/norm(r);  

if side=='N'
    %Gx=G*x
    
    Gx=zeros(4,length(x));
    iIn=reshape(iIn,4,[]); 
    
    for i2=1:size(ind,2)
        for i1=1:size(ind,1)
            
            for j0=1:4
                if iIn(j0,ind(i1,i2))
                    j1=i1+j1list(j0); j2=i2+j2list(j0); 
                    Gx(j0,ind(i1,i2))=r(1)*x(ind(i1,i2))+r(2)*x(ind(j1,j2));
                end
            end
            
        end
    end
    
    Gx=Gx(:);
    
elseif side=='T'
    %Gx=G'*x
    
    Gx=zeros(numel(ind),1);  
    iIn=reshape(iIn,4,numel(ind)); 
    x=reshape(x,4,numel(ind)); 
    
    for i2=1:size(ind,2)
        for i1=1:size(ind,1)
            for j0=1:4
                if iIn(j0,ind(i1,i2))
                    j1=i1+j1list(j0); j2=i2+j2list(j0); 
                    Gx(ind(i1,i2))=Gx(ind(i1,i2))+r(1)*x(j0,ind(i1,i2));
                    Gx(ind(j1,j2))=Gx(ind(j1,j2))+r(2)*x(j0,ind(i1,i2));
                end
            end
        end
        
    end
    
    Gx=Gx(:); 
    
end

end



%% 

function x=projectionG(ind,x,b,iIn)
%x=projectionG(ind,x,b,iIn). Calculates the projection of b on the space
%spanned by the rows iIn of G. 

i0=0; %iteration counter
J=[]; 
while true
    r=b-G(ind,'T',x,iIn); %calculate residual
    Gr=G(ind,'N',r,iIn);  %Calculate G'*r showing the orthogonality of the residual to row space G
    if i0==0; Gb=Gr; end

    %break conditions
    J=[J,max(Gr(iIn))]; 
    if i0>=20e3
        break;
    elseif norm(Gr(iIn))/norm(Gb(iIn))<1e-3
        break; 
    else
        if mod(i0,50)==0; display(sprintf('projectionG. Iteration %d. Max projection: %e',i0,J(end))); end
    end

    dir=Gr; %direction of greatest descent 
    Gd=G(ind,'T',dir,iIn); 
    
    gamma=norm(Gr)^2/norm(Gd)^2; %minimalization cost-function along dir 
    x=x+gamma*dir; %new estimate projection
    i0=i0+1;
end

end

