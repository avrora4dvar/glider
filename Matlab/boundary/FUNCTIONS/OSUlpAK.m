function [uf,tf]=OSUlpAK(u,t);
%clear all
%load tmp;
%u=wind_ustr;v=wind_vstr;t=time;
%
% [uf,tf]=OSUlpAK(u,t) ...
%			 is a 40 hour lowpass filter documented in
% the CODE-2 data report (gray cover, 8.5"x14", spiral bound, WHOI
% Technical Report 85-35), on page 21,
% The weights are of the form:
%   0.5*(1.+cos(pi*T/61))* (sin(2*pi*T/40))/(2*pi*T/40),
% where T = -60,-59,-58,...,59,60  (121 weights in all).
% Normalize so sum of 121 weights = 1.00
%
% Data should be hourly or else a linear interpolation will be performed
%
% AK (12/3/03): u can be matrix, last dim of u corr to time 
% No interpolation to original time

%if nargin<2;
% help OSUlp
% return;
%end;
 
nn=size(u);
ndim=length(nn);
if size(u,ndim)~=length(t) 
 error('last dim of u is not eq. to length of t');
end  

% reshape u into 2D file:
nv=prod(nn(1:end-1));  % <- number of rows
nt=nn(end);
u1=reshape(u,[nv nt]);  % <- reshape u into 2D file, assign to u1 

% Interpolate to hourly data 
tf=t(1):1/24:t(length(t));
n=length(tf);
ui=interp1(t,u1',tf,'linear');
ui=ui';  % <- [nv x n]

uf=zeros(nv,n);

	   if n<120,uf=ui;return;end

% define weight function
 T=-60:1:60;
 w=zeros(1,121); w(61)=1;
 w(1:60)=0.5*(1.+cos(pi*T(1:60)/61)).* ...
         (sin(2*pi*T(1:60)/40))./(2*pi*T(1:60)/40);
 w(62:121)=0.5*(1.+cos(pi*T(62:121)/61)).* ...
           (sin(2*pi*T(62:121)/40))./(2*pi*T(62:121)/40);
 w=w/sum(w);
whos
% Lowpass first 60 hours
 for i=1:60;
  ww=w(60-i+2:121); ww=ww/sum(ww);
  uf(:,i)=sum(repmat(ww,[nv 1]).*ui(:,1:i+60),2); 
end;
% Lowpass body of the time series
for i=61:n-60;
  uf(:,i)=sum(repmat(w,[nv 1]).*ui(:,i-60:i+60),2); 
end;
% Lowpass last 59 hours
for i=n-59:n;
  ww=w(1:n-i+61); ww=ww/sum(ww);
  uf(:,i)=sum(repmat(ww,[nv 1]).*ui(:,i-60:n),2); 
end;

uf=reshape(uf,[nn(1:end-1) n]);



