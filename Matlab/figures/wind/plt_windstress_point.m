%% Plot v-wind


tLim=[2380,2435]+datenum('2005-01-01'); 
romsDir='/home/aruba/vol2/ipasmans/Exp/Exp26_free/'; 


%% Read model time

wind.t=[]; 
wind.u=[]; 
wind.v=[];
wind.dir=romsDir; 
ncstart=[224,223,1]; 
ncend=[224,223,Inf]; 
nccount=ncend-ncstart+[1,1,1]; 
wind.grd_ind=[224,223]; 

%% Read wind stress

con=dir(romsDir); 
for iFile=1:length(con)
    if isempty(regexp(con(iFile).name,'ocean_his','once')); continue; end
    fname=fullfile(romsDir,con(iFile).name); 
    display(fname); 
    
    t1=ncread(fname,'ocean_time',[1],[Inf])/24/3600+datenum('2005-01-01'); 
    t1=reshape(t1,1,[]); 
    u1=ncread(fname,'sustr',ncstart-[1,0,0],[2,1,Inf]); 
    u1=reshape(mean(u1,1),1,[]); 
    v1=ncread(fname,'svstr',ncstart-[0,1,0],[1,2,Inf]); 
    v1=reshape(mean(v1,2),1,[]); 
    
    wind.t=cat(2,wind.t,t1); 
    wind.u=cat(2,wind.u,u1); 
    wind.v=cat(2,wind.v,v1); 
end
    
%% Save


save('/home/server/student/homes/ipasmans/Data/article_glider/wind_stress.mat','wind'); 
