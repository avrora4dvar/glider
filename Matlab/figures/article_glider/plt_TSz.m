%% Plot z-profile of glider
%clear all; clc; %close all; 
addpath('~/Matlab/roms/seawater/'); 
clear all; clc; 

%If first run activate
if true
    plots=[]; 
end

inFile='/home/server/student/homes/ipasmans/Data/article_glider/TS_glider_ana.mat';
figName='/home/server/student/homes/ipasmans/Figures/article_glider/tmp'; 
%modelName={'Free','Combined','Glider Only','Glider Fixed','Ens'}; 
modelName={'Surface','Glider','Combined','Free'}; 

lim.t=[2392,2413]+datenum('2005-01-01'); 
lim.lon=[-Inf,Inf]; 
lim.lat=[-Inf,Inf]; 
lim.z=[-Inf,Inf]; 

mode='rms'; 
field='salt'; 
zBin=[0:4:200]'; 


%% Load

load(inFile); 
model=model([1,2,3,6]); 

if strcmpi(field,'temp'); unit='C'; end
if strcmpi(field,'salt'); unit='ppt'; end
if strcmpi(field,'rho'); unit='kgm^{-3}'; end

%% Create density fields

obs.rho=sw_pden(obs.salt,obs.temp,obs.p,0); 
for iMod=1:length(model)
    model(iMod).rho=sw_pden(model(iMod).salt,model(iMod).temp,obs.p,0); 
end
 

%% Calculate values

obs.lon=repmat(obs.lon,[1,size(obs.depth,2)]); 
obs.lat=repmat(obs.lat,[1,size(obs.depth,2)]); 
obs.t=repmat(obs.t,[1,size(obs.depth,2)]); 

valN=nan(length(zBin)-1,length(model)+1);
val=nan(length(zBin)-1,length(model)+1);

for iz=1:length(zBin)-1
    in=obs.t>=min(lim.t) & obs.t<=max(lim.t) & obs.lon>=min(lim.lon) & obs.lon<=max(lim.lon) &...
        obs.lat>=min(lim.lat) & obs.lat<=max(loc.lat) &...
        obs.depth>=zBin(iz) & obs.depth<zBin(iz+1);
    in=in & ~isnan(obs.(field)); 
    for iMod=1:length(model)
        in=in & ~isnan(model(iMod).(field)); 
    end
    in=reshape(in,[],1); 
     
    if ismember(mode,{'profile'})
        if strcmpi(field,'temp'); xlabelStr='Temperature [C]'; end
        if strcmpi(field,'salt'); xlabelStr='Salinity [ppt]'; end
        valN(iz,1)=sum(in); 
        val(iz,1)=mean(obs.(field)(in));
        for iMod=1:length(model)
            valN(iz,iMod+1)=sum(in);
            val(iz,iMod+1)=mean(model(iMod).(field)(in));
        end
    elseif ismember(mode,{'rms'})
        xlabelStr=sprintf('RMSE [%s]',unit); 
        valN(iz,1)=sum(in);
        val(iz,1)=sqrt(mean( (obs.(field)(in)-mean(obs.(field)(in))).^2 )); 
        for iMod=1:length(model)
            valN(iz,iMod+1)=sum(in); 
            val(iz,iMod+1)=sqrt(mean( (model(iMod).(field)(in)-obs.(field)(in)).^2 )); 
            val(iz,iMod+1)=val(iz,iMod+1);
        end
    elseif ismember(mode,{'std'})
        xlabelStr=sprintf('Error standard deviation [%s]',unit); 
        valN(iz,1)=sum(in);
        for iMod=1:length(model)
            valN(iz,iMod+1)=sum(in);
            val(iz,iMod+1)=std( model(iMod).(field)(in)-obs.(field)(in) ); 
        end
    elseif ismember(mode,{'bias'})
        xlabelStr=sprintf('Bias [%s]',unit); 
        valN(iz,1)=sum(in); 
        for iMod=1:length(model)
            valN(iz,iMod+1)=sum(in); 
            val(iz,iMod+1)=mean( model(iMod).(field)(in)-obs.(field)(in) ); 
        end
    elseif ismember(mode,{'corr'})
        xlabelStr='Correlation'; 
        valN(iz,1)=sum(in); 
        for iMod=1:length(model)
            valN(iz,iMod+1)=sum(in);
            if sum(in)==0; continue; end
            val(iz,iMod+1)=corr(...
                reshape(model(iMod).(field)(in),[],1),...
                reshape(obs.(field)(in),[],1)...
                ); 
            %val(end-2,:)=NaN*val(end-2,:); 
        end
    elseif ismember(mode,{'skill'})
        xlabelStr=sprintf('Skill'); 
        valN(iz,1)=sum(in); 
        for iMod=1:length(model)
            dobs=abs( obs.(field)(in)-mean(obs.(field)(in)) ); 
            dmod=abs( model(iMod).(field)(in)-mean(model(iMod).(field)(in)) ); 
            drms=abs( model(iMod).(field)(in)-obs.(field)(in) ); 
            drms=drms-mean(drms); 
            
            valN(iz,iMod+1)=sum(in);
            val(iz,iMod+1)=1-mean(drms.^2)/mean( (dobs+dmod).^2 );  
        end
    end
end

%
for iMod=1:size(valN,2)
    val(valN(:,iMod)<.2*nanmedian(valN(:,iMod)),iMod)=NaN; 
end


%% Figure

figure(); hold on; 
set(gcf,'color','w','units','inches','paperunits','inches'); 
set(gcf,'papersize',[8.5,11]*.5,'position',[0,0,8.5,11]*.5,'paperposition',[0,0,8.5,11]*.5); 

%defaults
fs=10; 
set(0,'DefaultLineLineWidth',1.1);
markerStyle={'o','x','s','d','p','^'}; 
markerStyle={'','','','','','',''}; 
cmap=lines(length(model)); 


%% Plot
zBin=.5*zBin(1:end-1)+.5*zBin(2:end); 


clear hP; hold on; 
for i2=1:length(model)
    
    
    hP(i2)=plot(val(:,i2+1),zBin,['-',markerStyle{i2}],...
        'color',cmap(i2,:),'markerfacecolor',cmap(i2,:),...
        'displayname',modelName{i2},'linewidth',2);
end
if strcmpi(mode,'profile')
    hP(length(hP)+1)=plot(val(:,1),zBin,...
        'color','k','displayname','glider','linewidth',2); 
end

set(gca,'box','on','layer','top'); 
set(gca,'xgrid','on','ygrid','on'); 
set(gca,'ylim',[0,200],'YDir','reverse');
xlabel(xlabelStr); 
ylabel('Depth [m]'); 
if ismember(mode,{'bias'})
    set(gca,'xlim',[-1,1]*max(abs(get(gca,'xlim')))); 
elseif ismember(mode,{'rms','std'})
    xlim=get(gca,'xlim'); xlim(1)=0; 
    xlim=[0,2.5];
    set(gca,'xlim',xlim); 
elseif ismember(mode,{'corr'})
    xlim=[-1,1]; 
    set(gca,'xlim',xlim); 
elseif ismember(mode,{'skill'})
    xlim=[0,1];
    set(gca,'xlim',xlim);
end
if strcmpi(mode,'temp'); title('Temperature'); end
if strcmpi(mode,'salt'); title('Salinity'); end

hLeg=legend(hP); 
set(hLeg,'location','west','visible','on'); 
title(sprintf('%s-%s',datestr(lim.t(1),'mm/dd'),datestr(lim.t(2),'mm/dd'))); 

%% Stats

in=zBin>25; 
for i2=1:size(val,2)
   stat.rms(i2)=nansum(valN(in,i2).*val(in,i2).^2)/nansum(valN(in,i2));  
   %stat.rms(i2)=nanmean(val(in,i2).^2); 
end
stat.rms=sqrt(stat.rms)


%% save

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 18, 'fontWeight', 'bold');
%print(gcf,[figName,'_',mode,'_',field],'-dpng','-r200'); 

plot1=struct('mode',mode,'field',field,'lim',lim,'z',zBin,'val',val,'figName',figName,'name',{modelName});  
plots=[plots,plot1]; 
save('~/Data/article_glider/zprofile_plots_ensfor.mat','plots'); 








