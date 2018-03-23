function [ hs ] = smooth_contour(lon,lat,h,levels,dstep)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

close all; 
 
dl=5;  

x=(lon-mean(lon(:))).*cos(deg2rad(lat)); 
y=(lat-mean(lat(:))); 
x=deg2rad(x)*earthRadius;
y=deg2rad(y)*earthRadius;

[c1,c2]=contour(x,y,h,levels); 

V=[]; 
child=get(c2,'children');
for k=1:length(child)
    display(sprintf('Contour %d of %d',k,length(child))); 
    c1=get(child(k),'CData'); c1=c1(~isnan(c1)); 
    x1=get(child(k),'XData'); x1=x1(~isnan(x1)); 
    y1=get(child(k),'YData'); y1=y1(~isnan(y1)); 
    [unixy,xy2uni,uni2xy]=unique([x1,y1,c1],'rows','stable'); 
    x1=x1(xy2uni); y1=y1(xy2uni); c1=c1(xy2uni); 
    if length(x1)<2; continue; end
    
    
    d1=hypot(diff(x1),diff(y1)); 
    d1=[0;cumsum(d1)]; 
    d2=[0:dstep:d1(end)];
    if d2(end)<dstep*dl; continue; end
    
    x2=interp1(d1,x1,d2); 
    y2=interp1(d1,y1,d2); 
   
    
    if length(d2)<2*dl+1
        x3=mean(x2); y3=mean(y2); 
    elseif hypot(x1(end)-x1(1),y1(end)-y1(1))<dstep
        x3=nan(size(d2)); y3=nan(size(d2)); 
        for i1=1:length(x3)
            [p1,s,mu]=polyfit([-dl:dl],x2(1:2*dl+1),2); 
            p2=polyval(p1,[-dl:dl],s,mu); 
            x3(i1)=p2(1+dl); 
            [p1,s,mu]=polyfit([-dl:dl],y2(1:2*dl+1),2); 
            p2=polyval(p1,[-dl:dl],s,mu); 
            y3(i1)=p2(1+dl); 
            x2=circshift(x2,[0,-1]); y2=circshift(y2,[0,-1]); 
        end
    else
        x3=nan(1,length(d2)-2); y3=nan(1,length(d2)-2); 
        for i1=1:length(x3)
            bnd=[max(1,i1+1-dl):min(length(x2),i1+1+dl)]; 
            [p1,s,mu]=polyfit(bnd,x2(bnd),2); 
            p2=polyval(p1,bnd,s,mu); 
            x3(i1)=p2(bnd==i1+1); 
            [p1,s,mu]=polyfit(bnd,y2(bnd),2); 
            p2=polyval(p1,bnd,s,mu); 
            y3(i1)=p2(bnd==i1+1); 
        end
    end
        
    x3=reshape(x3,[],1); y3=reshape(y3,[],1); c3=c1(1)*ones(size(x3)); 
    V=[V;x3,y3,c3]; 
    
end

figure(); 
scatter(V(:,1),V(:,2),4*ones(size(V,1),1),V(:,3)); 

F=scatteredInterpolant(V(:,1),V(:,2),V(:,3)); 
hs=F(x,y); 


end

