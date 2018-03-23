function [sOut] = reconsruct_pcg(fileIn,flag_pre)
%RECONSTRUCT_PCG sOut=reconstruct_pcg(fileIn)
%   Reconstruct PCG output

info=ncinfo(fileIn); 


%% read observations

pcg.sig=ncread(fileIn,'sig'); 
pcg.obs=ncread(fileIn,'obs'); 
pcg.type=ncread(fileIn,'type');
typeUni=unique(pcg.type); 
    
precon=ones(size(type));
if flag_pre
    in=type<4; precon(in)=sum(in);
    in=type==4; precon(in)=sum(in);
    in=type>4; precon(in)=sum(in);
    precon=double(precon).^(-.25); 
end
if strcmpi('preDia',{info.Variables(:).Name})
    precon=ncread(fileIn,'preDia'); 
end



error('ok'); 


%% read residuals

r=ncread(fileIn,'r'); 
Br=ncread(fileIn,'Br');

n_iter=ncread(fileIn,'n_iter');
n_par=size(r,2)/(n_iter+1); 
sOut.n_iter=n_iter; sOut.n_par=n_par; 

rRes=r(:,end-n_par+1:end); 
r=r(:,1:end-n_par); Br=Br(:,1:end-n_par); 
Ar=Br+bsxfun(@times,r,precon.^2); 

r0=ncread(fileIn,'r0'); 
r00=bsxfun(@times,1./precon,r0); 

%% calculate T

rBr=r'*Br;
sOut.T=rBr\(Br'*Ar); 
sOut.rBr=rBr; 

%Calculate relative B-error
for k=1:n_par:size(rBr,1)
    sOut.Bnorm( floor(k/n_par)+1 )=sqrt(...
        sum( diag(rBr(k:k+n_par-1,k:k+n_par-1)) )/sum( diag(rBr(1:n_par,1:n_par)) )); 
end

%% Go through iterations


x=0*r0; Bx=0*r0; 
sOut.J(1,:)=.5*sum(x.*Bx,1)+.5*sum((r00-Bx).^2,1);
sOut.DAcor=zeros(size(obs,1),n_iter+0); 

for iIter=0:n_iter
   
    i2=n_par*iIter;
    
    if i2==0
        x=0*r0;
        Bx=0*r0;
    else
        rBr=r(:,1:i2)'*Br(:,1:i2);
        T=rBr\(Br(:,1:i2)'*Ar(:,1:i2));
        
        Pr0=rBr\(Br(:,1:i2)'*r0);
        Px=T\Pr0;
        x=r(:,1:i2)*Px; x=bsxfun(@times,x,precon); 
        Bx=Br(:,1:i2)*Px; Bx=bsxfun(@times,Bx,1./precon); 
        
        sOut.DAcor(:,iIter+1)=sig.*Bx(:,1); 
    end
    
    %Change in x in B-norm
    if iIter==0
        sOut.diffx(1,:)=nan(1,size(x,2)); 
    else
        sOut.diffx(iIter+1,:)=sqrt( sum((x-xprior).*(Bx-Bxprior),1)./sum(x.*Bx,1) ); 
    end
    xprior=x; Bxprior=Bx; 
      
    %Cost function 
    sOut.J(iIter+1,:)=.5*sum(x.*Bx,1)+.5*sum((r00-Bx).^2,1);
    
    %rms
    for it=1:length(typeUni)
        in=type==typeUni(it); 
        sOut.iterRms.(sprintf('type%-.4d',typeUni(it)))(iIter+1,:)=rms(r00(in,:)-Bx(in,:),1); 
    end
    sOut.iterRms.all(iIter+1,:)=rms(r00-Bx,1); 
    
end


end

