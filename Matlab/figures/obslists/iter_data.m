function [stat,iter,obs] = iter_data(fname,varargin)
%PCG_DATA Summary of this function goes here
%   Detailed explanation goes her

if isempty(varargin)
    iMem=1; nMem=Inf; 
else
    iMem=varargin{1}; nMem=1; 
end

%% Read data 

obs.sig=ncread(fname,'sig'); 
obs.type=ncread(fname,'type'); 
obs.val=ncread(fname,'obs'); 
obs.valn=obs.val./obs.sig; 

%Search directions
r=ncread(fname,'r'); 
Br=ncread(fname,'Br'); 

%Input
r0=ncread(fname,'r0',[1,iMem],[Inf,nMem]); 
obs.for=bsxfun(@times,r0,obs.sig); 
obs.for=bsxfun(@minus,obs.val,obs.for); 
% 
% %Preconditoner
% pre.U=ncread(fname,'preU'); 
obs.preS=ncread(fname,'preS'); 

%Number of iterations
nMax=ncread(fname,'n_iter'); %iterations 0,1,...,nmax

%Number parallel
nPar=size(r,2)/(nMax+1); 

r=r(:,1:nPar*nMax); Br=Br(:,1:nPar*nMax); 
rBr=r'*Br; 

%% Reconstruct each iteration
 
for n=0:nMax
   display(sprintf('Processing iteration %d of %d',n,nMax)); 
   iter1=struct('J',[],'x',[],'Bx',[],'rBr',[],'T',[],'range',[],'nPar',nPar,'n',n); 
    
   if n==0
       iter1.range=[0,0]; 
       iter1.x=0*r0; iter1.Bx=0*r0; 
   else
       %Search vectors
       iter1.range=[1:n*nPar]; 
       r1=r(:,iter1.range); Br1=Br(:,iter1.range); 
       Ar1=Br1+r1; 
       
       %Projections
       iter1.rBr=r1'*Br1;
       Pr=iter1.rBr\(Br1'*r0);
       iter1.T=iter1.rBr\(Br1'*Ar1);
       
       %Find solution
       iter1.x=r1*(iter1.T\Pr);
       iter1.Bx=Br1*(iter1.T\Pr);
       
       %check x
       xAlt=(Br1'*Ar1)\(Br1'*r0); xAlt=r1*xAlt; 
       for i2=1:size(xAlt,2)
           if norm(iter1.x(:,i2)-xAlt(:,i2))/norm(xAlt(:,i2))>1e-3
               display('error x'); 
           end
       end
      
   end
   
   %Find Ritz values
   if ~isempty(iter1.T)
       [iter1.lanczos_vector,S]=eig(iter1.T); 
       iter1.ritz=diag(S); 
       if max(abs(imag(iter1.ritz)))>1e-3*min(abs(iter1.ritz))
           display('Imaginary Ritz values'); 
           iter1.ritz=iter1.ritz; 
       else
           iter1.ritz=real(iter1.ritz); 
       end
   end
       
   %Find normalized residual
   Ax=iter1.Bx+iter1.x; 
   iter1.r=r0-Ax; 
   
   %Find analysis
   iter1.ana=bsxfun(@times,obs.sig,iter1.Bx); 
   iter1.ana=obs.for+iter1.ana; 
      
   %Find cost-function
   for i2=1:size(r0,2)
      iter1.J(i2)=.5*(iter1.x(:,i2)'*iter1.Bx(:,i2))+...
          .5*norm(r0(:,i2)-iter1.Bx(:,i2))^2; 
   end
   
   %Find norm residual
   P1=rBr\(Br'*iter1.r);
   for i2=1:size(r0,2)
      iter1.rnorm(i2)=norm(iter1.r(:,i2).*obs.sig); 
      iter1.rnormN(i2)=norm(iter1.r(:,i2)); 
      iter1.rnormB(i2)=sqrt(P1(:,i2)'*rBr*P1(:,i2)); 
   end
   
   %Change x with respect to previous iteration
   if n==0
       iter1.xnorm=nan(1,size(r0,2)); 
       iter1.xnormN=nan(1,size(r0,2)); 
       iter1.xnormB=nan(1,size(r0,2)); 
   else
      dx=iter1.x-iter{n}.x;
      P1=rBr\(Br'*dx); P2=rBr\(Br'*iter1.x); 
      for i2=1:size(r0,2)
         iter1.xnorm(i2)=norm(dx(:,i2).*obs.sig)/norm(iter1.x(:,i2).*obs.sig); 
         iter1.xnormN(i2)=norm(dx(:,i2))/norm(iter1.x(:,i2));
         iter1.xnormB(i2)=sqrt(P1(:,i2)'*rBr*P1(:,i2))/sqrt(P2(:,i2)'*rBr*P2(:,i2)); 
      end
   end
   
   %Orthogonality unobserved/observed part
   for i2=1:size(r0,2)
       iter1.ortho(i2)=(iter1.x(:,i2)'*iter1.Bx(:,i2))/norm(iter1.x(:,i2))/norm(iter1.Bx(:,i2)); 
   end
   iter1.ortho=rad2deg(acos(iter1.ortho)); 

   %Set output
   iter{n+1}=iter1; 
    
end

%% Collect stats

display('Collecting from iterations'); 

for n=1:length(iter)
    stat.J(n,:)=iter{n}.J; 
    stat.rnorm.sig(n,:)=iter{n}.rnorm./iter{1}.rnorm; 
    stat.rnorm.N(n,:)=iter{n}.rnormN./iter{1}.rnormN; 
    stat.rnorm.B(n,:)=iter{n}.rnormB./iter{1}.rnormB; 
    stat.xnorm.sig(n,:)=iter{n}.xnorm; 
    stat.xnorm.N(n,:)=iter{n}.xnormN; 
    stat.xnorm.B(n,:)=iter{n}.xnormB; 
    stat.xnorm.P(n,1)=2*(sqrt(max(iter{end}.ritz)/min(iter{end}.ritz))-1)^(n)/...
        (sqrt(max(iter{end}.ritz)/min(iter{end}.ritz))+1)^(n); 
    stat.ortho(n,:)=iter{n}.ortho; 
    
    %Norm per type
    typeUni=unique(obs.type); 
    for i1=1:length(typeUni)
        in=obs.type==typeUni(i1); 
        for i2=1:size(r0,2)
        stat.rms.(sprintf('type%-.4d',typeUni(i1)))(n,i2)=rms(iter{n}.r(in,i2)); 
        end
        stat.rms.all(n,i2)=rms(iter{n}.r(:,i2)); 
    end
    
    if n>1
        stat.ritz.max(n)=max(iter{n}.ritz);
        stat.ritz.min(n)=min(iter{n}.ritz);
    else
        stat.ritz.max(n)=NaN; stat.ritz.min(n)=NaN;
    end
end
end


