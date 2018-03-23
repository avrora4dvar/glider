function [z] = roms_sigma2z(sigma,h,varargin)
%roms_sigma2z converst sigma to z-coordinates
%  roms_sigma2z(sigma,h,zeta,par) converts sigma-coordinates sigma (1D array) to 
%  z-coordinates using bottom depth h (2D array_, sea surface height zeta (2D aray) and grid
%  parameters par. The latter 2 arguments are optional. If no zeta is
%  specified zeta=0 is assumed. 

%% Input

%Grid parameter
if ~isempty(varargin) && isstruct(varargin{end})
    par=varargin{end};
else
    par.Vtransform=2;
    par.Vstretching=4;
    par.Tcline=50;
    par.theta_s=8;
    par.theta_b=3; 
    
    par.alpha=NaN; par.beta=NaN; par.gamma=NaN; 
end

%zeta
if ~isempty(varargin) && ~isstruct(varargin{1})
    zeta=varargin{1};
else
    zeta=zeros(size(h)); 
end

%% Initiate output

%sigma coordinates
sigma=repmat(reshape(sigma,1,1,[]),[size(h,1),size(h,2),1]); 

%bottom/surface
h=repmat(h,[1,1,size(sigma,3)]); 
zeta=repmat(zeta,[1,1,size(sigma,3)]); 

%% Calculate

%hc
if par.Vtransform==1
    hc=min(min(h(:)),par.Tcline); 
else
    hc=par.Tcline;
end

%stretching function
C=roms_Vstretching(par.Vstretching,sigma,par.theta_s,par.theta_b,par.alpha,par.beta,par.gamma); 

%output
z=roms_Vtransform(par.Vtransform,sigma,C,h,zeta,hc); 

end

%% Vtransform

function [z] =roms_Vtransform(type,sigma,C,h,zeta,hc)
%ROMS_VTRANSFORM roms_Vtransform calculates z coordinates for sigma-grid
%   roms_Vtransform calculates z-coordinates given the value of the
%   parameter Vtransform (type), the sigma coordinates (sigma), the stretching function of the layers (C), bottom
%   depth h, and sea surface height w.r.t. reference (zeta)

switch type
    case 1 %Vtransform==1
        S=hc*sigma+(h-hc).*C;
        z=S+zeta.*(1+S./h); 
    case 2
        S=(hc*sigma+h.*C)./(hc+h); 
        z=zeta+(zeta+h).*S; 
    otherwise
        error('argument type must be either 1 or 2'); 
end


end

%% Vstretching

function C=roms_Vstretching(type,sigma,theta_s,theta_b,alpha,beta,gamma)
%ROMS_VSTRETCHING

switch type
    case 1 %song and haidvogel (1994)
        C=(1-theta_b)*sinh(theta_s*sigma)./sinh(theta_s)+.5*theta_b*(tanh(theta_s*(sigma+.5))./tanh(.5*theta_s)-1); 
    case 2 %Shchepetkin (2005)
        mu=(sigma+1).^alpha*(1+alpha/beta*(1-(sigma+1).^beta)); 
        Csur=(1-cosh(theta_s*sigma))/(cosh(theta_s)-1);
        Cbot=sinh(theta_b*(sigma+1))/sinh(theta_b)-1; 
        C=mu.*Csur+(1-mu).*Cbot; 
    case 3 %Geyer
        mu=.5*(1-tanh(gamma*(sigma+.5))); 
        Csur=-log(cosh(gamma*abs(sigma).^theta_s))/log(cosh(gamma)); 
        Cbot=log(cosh(gamma*(sigma+1).^theta_b))/log(cosh(gamma))-1;
        C=mu.*Cbot+(1-mu).*Csur; 
    case 4 %Shchepetkin (2010)
        if theta_s>0
            C=(1-cosh(theta_s*sigma))/(cosh(theta_s)-1); 
        else
            C=-sigma.^2;
        end
        C=(exp(theta_b*C)-1)/(1-exp(-theta_b)); 
    otherwise
        error('type must be 1,2,3 or 4'); 
end
end