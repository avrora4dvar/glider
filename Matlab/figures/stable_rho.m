function [rho] = stable_rho(rho,dz,varargin)
%STABLE_RHO [rho]=stable_rho(rho,dz)
%   Transforms density profile into stable density profile whilst
%   conserving mass. rho,dz are mxnxk-matrices

s=size(rho,1)*size(rho,2); 

%Minimal density difference between different layers
iPar=find(strcmpi('min_step',varargin),1); 
if isempty(iPar)
    dmin=0; 
else
    dmin=varargin{iPar+1}; 
end

%Find unstable profiles
drho=diff(rho,1,3);
mask=-max(drho,[],3)>dmin | any(isnan(drho),3);

%Create stable profile
for i1=1:size(rho,1)
    for i2=1:size(rho,2)
        if mask(i1,i2); continue; end
        
        rho(i1,i2,:)=stable_rho1(squeeze(rho(i1,i2,:)),squeeze(dz(i1,i2,:)),dmin);
        mask(i1,i2)=true;
    end
end
        
%% Stabilization full profile

% for i1=size(rho,3):-1:2 %shallow layer
%     for i2=i1-1:-1:1 %deeper layer
%         
%         %Density difference between deep layer and shallow layer must be
%         %sufficiently positive
%         drho=diff(rho(:,:,[i1,i2]),1,3); 
%         in=find(drho<dmin*abs(i1-i2)); 
%         if isempty(in); continue; end
%         
%         %Calculate average density over the layers
%         rho1=sum(rho(:,:,[i2:i1]).*dz(:,:,[i2:i1]),3)-...
%             sum(bsxfun(@times,dz(:,:,[i2:i1]),-dmin*reshape([i2:i1]-i1,1,1,[])),3); 
%         rho1=rho1./sum(dz(:,:,[i2:i1]),3);
%         
%         %Set new stable density
%         for i3=i2:i1
%             rho(in+(i3-1)*s)=rho1(in)+dmin*abs(i3-i1); 
%         end
%     end
% end


end

%% One dimensional function

function rho=stable_rho1(rho,dz,dmin)
    
    for i1=length(dz):-1:2
        for i2=i1-1:-1:1
            %Difference
            drho=rho(i2)-rho(i1); 
            if drho>dmin*(i1-i2); continue; end
            %Calculate new average density
            rho1=sum(rho(i2:i1).*dz(i2:i1))+...
                dmin*sum(dz(i2:i1).*([i2:i1]'-i1)); 
            rho1=rho1/sum(dz(i2:i1)); 
            %Set new density
            rho(i2:i1)=rho1-([i2:i1]-i1)*dmin; 
        end
    end

end

