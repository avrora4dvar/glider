function [in0,iclosest]=find_closest(XX,YY,mask)

in0=find(mask==0); % 1d arrays
in1=find(mask==1);
n0=length(in0);
if n0>0
 itmp=zeros(n0,1);
 for k=1:n0
  [tmp,itmp(k)]=min( (XX(in1)-XX(in0(k))).^2+(YY(in1)-YY(in0(k))).^2);
 end
 iclosest=in1(itmp);
else
 iclosest=[];
end % n0>0
