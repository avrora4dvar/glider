function [roms_index,t]=get_hash_roms(roms_dir,tLim,varargin)
%GETHASHROMS Get list of ocean_his numbers and the times stored in them

dateRef=datenum('2005-01-01'); 
hashTable=[]; 

if isempty(varargin)
    lineBase='ocean_his'; 
else
    lineBase=varargin{1}; 
end

if ~exist(roms_dir,'dir')
    error('%s does not exist',roms_dir); e
end

dirCon=dir(roms_dir); 
for i0=1:length(dirCon)
   line=regexp(dirCon(i0).name,[lineBase,'_(\d+).nc'],'tokens'); 
   if ~isempty(line)
      line=line{1}; 
      fileNo=str2num(line{1}); 
      
      %size oceantimec
      info=ncinfo(fullfile(roms_dir,dirCon(i0).name),'ocean_time'); 
      dim_t=info.Size;
      
      %times in file
      tBnd(1)=ncread(fullfile(roms_dir,dirCon(i0).name),'ocean_time',[1],[1]); 
      tBnd(2)=ncread(fullfile(roms_dir,dirCon(i0).name),'ocean_time',[dim_t],[1]); 
      t=linspace(tBnd(1),tBnd(2),dim_t)/24/3600+dateRef; 
      
      %write to table
      tmp=[ones(length(t),1)*fileNo,[1:dim_t]',t(:)]; 
      hashTable=[hashTable;tmp]; 
   end
end


%select unique times
tLim(1)=max(min(hashTable(:,3)),tLim(1)); 
tLim(2)=min(max(hashTable(:,3)),tLim(2));

iTable(1)=find(hashTable(:,3)<=tLim(1),1,'last'); 
iTable(2)=find(hashTable(:,3)>=tLim(2),1,'first'); 
hashTable=hashTable(iTable(1):iTable(2),:); 

hashTable=sortrows(hashTable,[3,2]); %Takes for the analysis the first field of the window
[tUni,t2uni,uni2t]=unique(hashTable(:,3)); 
hashTable=hashTable(t2uni,:);
[fUni,f2uni,uni2f]=unique(hashTable(:,1)); 
for i0=1:length(fUni)
    roms_index(i0).filename=fullfile(roms_dir,sprintf([lineBase,'_%-.4d.nc'],fUni(i0))); 
    roms_index(i0).index=hashTable(hashTable(:,1)==fUni(i0),2); 
    roms_index(i0).time=hashTable(hashTable(:,1)==fUni(i0),3); 
end

%number of times
t=tUni; 

end

