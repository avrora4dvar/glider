%% Calculate enstrophy
clear all; clc; 

inFile={...
    '/home/server/student/homes/ipasmans/Data/article_glider/Exp35ana_his.mat',...
    '/home/server/student/homes/ipasmans/Data/article_glider/Exp36ana_his.mat',...
    '/home/server/student/homes/ipasmans/Data/article_glider/Exp37ana_his.mat',...
    '/home/server/student/homes/ipasmans/Data/article_glider/Exp40ana_his.mat'...
    }; 

% inFile={...
%     '/home/server/student/homes/ipasmans/Data/article_glider/Exp40ana_his.mat',...
%     '/home/server/student/homes/ipasmans/Data/article_glider/Exp37ana_his.mat',...
%     '/home/server/student/homes/ipasmans/Data/article_glider/Exp36ana_his.mat',...
%     '/home/server/student/homes/ipasmans/Data/article_glider/Exp43ana_his.mat',...
%     '/home/server/student/homes/ipasmans/Data/article_glider/Exp36aana_his.mat'
%     }; 

modName={'Exp35ana','Exp36ana','Exp37ana','Exp40'}; 
flag_plot=false; 
xtick=[2392:3:2413]+datenum('2005-01-01'); 

%% Read grid

grdFile='/home/aruba/vol2/ipasmans/Exp/Prm/grd_ow2km_r13lin_mix.nc';

grd.lon=ncread(grdFile,'lon_rho'); 
grd.lat=ncread(grdFile,'lat_rho'); 
grd.pn=ncread(grdFile,'pn'); 
grd.pm=ncread(grdFile,'pm'); 
grd.A=1./(grd.pn.*grd.pm); 

grd.dz=diff( ncread(grdFile,'z0_w'),[],3 );

%Select region
grd.mask=zeros(size(grd.lon));
grd.mask(grd.lon>-127 & grd.lat<47 & grd.lat>42)=1; 
grd.A(~grd.mask)=0; 

%% 

for iMod=1:length(inFile)
    mod1=load(inFile{iMod}); 
    
    %Calculate enstrophy
    enstr=.5*bsxfun(@times,grd.A(2:end-1,2:end-1),mod1.vort.^2);
    model{iMod}.enstrophy=nansum(reshape(enstr,[],size(enstr,3)),1); 
    model{iMod}.t=mod1.t; 
    
    %Fit rate
    for it=1:length(xtick)-1
        in=model{iMod}.t>xtick(it) & model{iMod}.t<=xtick(it+1); 
        [B,BINT,R,RINT,STATS]=regress(model{iMod}.enstrophy(in)',...
                [ones(sum(in),1),(model{iMod}.t(in)-xtick(it))'*24*3600]);
        model{iMod}.slope.t(it)=mean(xtick(it:it+1));  
        model{iMod}.slope.M(it)=B(2); 
        model{iMod}.slope.U(it)=B(2)+BINT(2,2); 
        model{iMod}.slope.L(it)=B(2)-BINT(2,1); 
    end
        
    
    clear mod1; 
end

%% Save

if ~flag_plot
   save('/home/server/student/homes/ipasmans/Data/article_glider/enstrophy_ana.mat','model','inFile'); 
   error('ok')
end

%% Plot

if flag_plot
    
    close all; figure(); hold on;
    cmap=lines(length(model));
    for iMod=1:length(model)
        hP(iMod)=plot(model{iMod}.t,model{iMod}.enstrophy,'o-','color',cmap(iMod,:),'displayname',modName{iMod},'linewidth',1.5);
    end
    
    xtick=[datenum('2011-07-21 00:00'):3:datenum('2011-08-11 00:00')];
    set(gca,'xgrid','on','ygrid','on','box','on');
    set(gca,'xtick',xtick,'xticklabel',cellstr(datestr(xtick,'mm/dd')));
    xlabel('2011'); ylabel('Surface layer enstrophy [m^2 s^{-2}]');
    
    legend(hP,'location','northwest');
    print(gcf,'~/Figures/article_glider/enstrophy','-dpng','-r300');
    
end

%% Slope

if flag_plot
    
    slope=nan(length(model),length(xtick)-1);
    for it=1:length(xtick)-1
        for iMod=1:length(model)
            in=model{iMod}.t>xtick(it) & model{iMod}.t<xtick(it+1);
            p1=polyfit((model{iMod}.t(in)-xtick(it))*24*3600,model{iMod}.enstrophy(in),1);
            
            
            [B,BINT,R,RINT,STATS]=regress(model{iMod}.enstrophy(in)',...
                [ones(sum(in),1),(model{iMod}.t(in)-xtick(it))'*24*3600]);
            
            slope(iMod,it)=B(2);
            slopeL(iMod,it)=BINT(2,1);
            slopeU(iMod,it)=BINT(2,2);
        end
    end
    
    figure(); hold on;
    for iMod=1:length(model)
        hP(iMod)=errorbar(.5*xtick(1:end-1)+.5*xtick(2:end),...
            slope(iMod,:),...
            slope(iMod,:)-slopeL(iMod,:),slopeU(iMod,:)-slope(iMod,:),...
            'o-','color',cmap(iMod,:),'displayname',modName{iMod},'linewidth',1.5);
    end
    
    set(gca,'xgrid','on','ygrid','on','box','on');
    set(gca,'xtick',xtick,'xticklabel',cellstr(datestr(xtick,'mm/dd')));
    xlabel('2011'); ylabel('Rate of change enstrophy [m^2 s^{-3}]');
    legend(hP,'location','northwest');
    print(gcf,'~/Figures/article_glider/enstrophy_rate','-dpng','-r300');
    
end

%% Stat test


%Collect values
val=[]; 
for iMod=1:length(model)
    in=model{iMod}.t>=min(xtick) & model{iMod}.t<=max(xtick); 
    val=[val; model{iMod}.enstrophy(in)]; 
    valt=model{iMod}.t(in); 
end

val=(val(:,2:end)-val(:,1:end-1)).^2; 
in=[3:3:length(valt)-1]; 

%Between window
jump=sqrt(sum(val(:,in),2)/sum(in));
%Withtin window
njump=sqrt(sum(val(:,~in),2)/sum(~in)); 
