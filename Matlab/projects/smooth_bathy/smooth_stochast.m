function [hs,ch,classes] = smooth_stochast(h,rmax)
%smooth_max Summary of this function goes here
%   Detailed explanation goes here

%NaN values
hnan=3; 
h(isnan(h))=hnan; h=max(h,hnan); 
h0=h;
T=[1,1]; 

%convert to log
H=h; 
h0=log(h0/hnan); 
h=log(h/hnan); 


%classify
classes=[0:.5:99]*rmax; 
classes=classes(1:find(classes>max(h(:)),1,'first'));
ch=zeros(size(h)); ch0=zeros(size(h)); 
for k=1:length(classes)-1
    in=h0>classes(k) & h0<=classes(k+1); 
    ch0(in)=k; 
    ch(in)=k; 
end

%samples

for k=1:1e6
    if mod(k,1e4)==0
        display(sprintf('k=%d',k))
    end
    
    i1=interp1([0,1],[1.51,size(ch,1)-.51],rand()); i1=round(i1); 
    i2=interp1([0,1],[1.51,size(ch,2)-.51],rand()); i2=round(i2); 
    
%     %Select range possible new classes
    cN=[ch(i1-1,i2),ch(i1+1,i2),ch(i1,i2-1),ch(i1,i2+1),...
        ch(i1-1,i2-1),ch(i1-1,i2+1),ch(i1+1,i2-1),ch(i1+1,i2+1)]; %possible classes
    hN=[H(i1-1,i2),H(i1+1,i2),H(i1,i2-1),H(i1,i2+1),...
        H(i1-1,i2-1),H(i1-1,i2+1),H(i1+1,i2-1),H(i1+1,i2+1)]; 
    cUni=unique(cN(1:4)); 
   
    
    %Draw proposal
    j1=floor(rand()*length(cUni))+1;
    c1=cUni(j1); 
    
    %Prob from old to new state
    g0=1/length(cUni); 
    %Prob from new to old stat
    g1=any(cN==ch(i1,i2))/length(cUni); 
    
    %Energy different class
    dE1=abs(c1-ch0(i1,i2)).^2; 
             
    %Energy different surface
    dE2=surface_energy(c1,cN)-surface_energy(ch(i1,i2),cN); 
    
    
    alpha=exp(-dE1/T(1)-dE2/T(2))*g1/g0;
    [alpha,exp(-dE1/T(1)-dE2/T(2)),g0,g1];
    if alpha>rand()
        H(i1,i2)=mean(hN(cN==c1)); 
        ch(i1,i2)=c1; 
    end
    
    
    
end



%convert back to depth
hs=H; 
hs(hs<=hnan)=NaN; 



end

function E=surface_energy(c0,c1)
    E=0; 
    for k=1:length(c1)
        if c1(k)==c0
            E=E; 
        elseif c1(k)==0 || c0==0
            E=E+1; 
        else
            E=E+exp(max(0,abs(c1(k)-c0)-1)); 
        end
    end
end



