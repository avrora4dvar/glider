% find daily averages for years 1969-2013 (over 24years)
% for Columbia River DISCHARGE
rd=0;
if rd==1,
fname='CR/Discharge_1969_2013.txt';
fid=fopen(fname,'r');
for k=1:28,s=fgetl(fid);end
date=[];D=[];
while length(s)>0,
s=fgetl(fid);
if s==-1, break;end
if s(1:4)~='USGS',break;end
ik=findstr(s,'-');
if isempty(ik)>0,break;end
id=findstr(s,'A');
if isempty(id)==0,
 date=[date;s(ik(1)-4:ik(1)+5)];
 D=[D,str2num(s(ik(2)+3:id(1)-1))];
end
end
fclose(fid);
end
% make daily averages over years
yyyy=str2num(date(:,1:4));mm=str2num(date(:,6:7));dd=str2num(date(:,9:10));
% for a NON-visocos year, say 2001
% For Feb 29 will take ave btw March 1 and Feb 28
t0=datenum(2001,1,1);
Da=zeros(size(365));La=Da;
for iday=1:365
  d1=datestr(t0+iday-1,'yyyy-mm-dd');
  y1=str2num(d1(:,1:4));mm1=str2num(d1(:,6:7));dd1=str2num(d1(:,9:10));
  ia=find(mm==mm1 & dd==dd1);
  Da(iday)=mean(D(ia));La(iday)=length(ia);
end
save CR_Dave.mat Da;
figure(1);clf
% Exclude Feb 29
inv=find(dd~=29 | mm~=02);
Di=D(inv);datei=date(inv,:);
yyyyi=yyyy(inv);ddi=dd(inv);mmi=mm(inv);
for y1=2010:2013
   d0=datenum(y1,1,1);
   iy=find(yyyyi==y1);
   D1=Di(iy);
   t1=datenum(yyyyi(iy),mmi(iy),ddi(iy))-d0;
   plot(t1,D1,'k');hold on
end
plot(Da,'Color','r','LineWidth',2);
