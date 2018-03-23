% find daily averages for years 1968-2003 (over 10-14years for various dates)
rd=0;
if rd==1,
fname='CR/Temp_1968_2002.txt';
fid=fopen(fname,'r');
for k=1:28,s=fgetl(fid);end
date=[];T=[];
while length(s)>0,
s=fgetl(fid);
if s==-1, break;end
if s(1:4)~='USGS',break;end
ik=findstr(s,'-');
if isempty(ik)>0,break;end
id=findstr(s,'A');
if isempty(id)==0,
 date=[date;s(ik(1)-4:ik(1)+5)];
 T=[T,str2num(s(ik(2)+3:id(1)-1))];
end
end
fclose(fid);
end
% make daily averages over years
yyyy=str2num(date(:,1:4));mm=str2num(date(:,6:7));dd=str2num(date(:,9:10));
% for a NON-visocos year, say 2001
% For Feb 29 will take ave btw March 1 and Feb 28
t0=datenum(2001,1,1);
Ta=zeros(size(365));La=Ta;
for iday=1:365
  d1=datestr(t0+iday-1,'yyyy-mm-dd');
  y1=str2num(d1(:,1:4));mm1=str2num(d1(:,6:7));dd1=str2num(d1(:,9:10));
  ia=find(mm==mm1 & dd==dd1);
  Ta(iday)=mean(T(ia));La(iday)=length(ia);
end
save CR_Tave.mat Ta;

