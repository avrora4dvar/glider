%% read surface from files, LP-filter, interpolate it to rho-points if necessary and save to file with a greater time step

clear all; clc; close all; 

addpath('/home/server/student/homes/ipasmans/Matlab/roms/seawater/'); 
addpath('/home/server/student/homes/ipasmans/Matlab/roms/t_tide/'); 
addpath('..'); 

%% input

%grid
grid_file='/home/aruba/vol2/ipasmans/Exp/Prm/grd_ow2km_r13lin_mix.nc'; 
%layer to be read from ROMS files
layer=40; 
%times to save to output
t=[2392.5:2412.5]+datenum('2005-01-01'); 
%roms dir
roms_dir='/home/aruba/vol2/ipasmans/Exp/Exp35/Exp38_ana'
%roms_dir='/home/aruba/vol2/ipasmans/Exp/Exp29/Exp26a/'; 
roms_ref=[]; 
%output file
out_file={'/home/server/student/homes/ipasmans/Data/article_glider/Exp38ana_avg.mat',...
    ''}; 
%fields in output
fieldOut={'rho','temp','salt','u','v','vort'}; 
fieldOut={'u','v','vort','temp','salt'}; 
%filter
%filter='lanczos'; fcut=24/36; fwidth=48/24; 
%filter='step'; fwidth=12/24; 


%% read grid

%rho-grid
grd.rho.lon=ncread(grid_file,'lon_rho'); 
grd.rho.lat=ncread(grid_file,'lat_rho'); 
grd.rho.z=ncread(grid_file,'z0_r'); 
grd.rho.mask=ncread(grid_file,'mask_rho'); 

%u-grid
grd.u.lon=ncread(grid_file,'lon_u'); 
grd.u.lat=ncread(grid_file,'lat_u'); 
grd.u.mask=ncread(grid_file,'mask_u'); 
grd.u.latRad=repmat(deg2rad(grd.u.lat),[1,1,length(t)]); 
grd.u.lonRad=repmat(deg2rad(grd.u.lon),[1,1,length(t)]); 

%v-grid
grd.v.lon=ncread(grid_file,'lon_v'); 
grd.v.lat=ncread(grid_file,'lat_v'); 
grd.v.mask=ncread(grid_file,'mask_v'); 
grd.v.latRad=repmat(deg2rad(grd.v.lat),[1,1,length(t)]); 
grd.v.lonRad=repmat(deg2rad(grd.v.lon),[1,1,length(t)]); 

%p
grd.p.lat=.25*grd.rho.lat(1:end-1,1:end-1)...
    +.25*grd.rho.lat(2:end,1:end-1)...
    +.25*grd.rho.lat(1:end-1,2:end)...
    +.25*grd.rho.lat(2:end,2:end); 
grd.p.latRad=repmat(deg2rad(grd.p.lat),[1,1,length(t)]); 

s=size(grd.rho.mask)-[2,2]; 


%% read 2D ROMS output
clear val; clear val1; clear u; clear v;

fieldNames={'zeta','ubar','vbar','vortbar'}; 
for iField=1:length(fieldNames)
    if ~ismember(fieldNames{iField},fieldOut); continue; end
    
    display(sprintf('reading %s',fieldNames{iField})); 
    

    loc=struct('t',t,'layer',[]);
    if strcmpi(fieldNames{iField},'ubar')
        roms=get_roms_layer(roms_dir,fieldNames{iField},loc,'delta');
        if ~isempty(roms_ref)
            romsRef=get_roms_layer(roms_ref,fieldNames{iField},loc,'delta');         
            val1=roms.val-romsRef.val;
        else
            val1=roms.val; 
        end
        val=.5*val1(1:end-1,2:end-1,:)+.5*val1(2:end,2:end-1,:);
    elseif strcmpi(fieldNames{iField},'vbar')
        roms=get_roms_layer(roms_dir,fieldNames{iField},loc,'delta');
        if ~isempty(roms_ref)
            romsRef=get_roms_layer(roms_ref,fieldNames{iField},loc,'delta'); 
            val1=roms.val-romsRef.val;
        else
            val1=roms.val; 
        end
        val=.5*val1(2:end-1,1:end-1,:)+.5*val1(2:end-1,2:end,:);
    elseif strcmpi(fieldNames{iField},'vortbar')
        u=get_roms_layer(roms_dir,'ubar',loc,'delta');
        v=get_roms_layer(roms_dir,'vbar',loc,'delta');
        
        dv=(v.val(2:end,:,:)-v.val(1:end-1,:,:))./...
            (grd.v.lonRad(2:end,:,:)-grd.v.lonRad(1:end-1,:,:)); 
        du=(u.val(:,2:end,:).*cos(grd.u.latRad(:,2:end,:))...
            -u.val(:,1:end-1,:).*cos(grd.u.latRad(:,1:end-1,:)))...
            ./(grd.u.latRad(:,2:end,:)-grd.u.latRad(:,1:end-1,:)); 
        val1=(dv-du)./grd.p.latRad/earthRadius; 
      
       
        if ~isempty(roms_ref)
            u=get_roms_layer(roms_ref,'ubar',loc,'delta');
            v=get_roms_layer(roms_ref,'vbar',loc,'delta');
            
            dv=(v.val(2:end,:,:)-v.val(1:end-1,:,:))./...
                (grd.v.lonRad(2:end,:,:)-grd.v.lonRad(1:end-1,:,:));
            du=(u.val(:,2:end,:).*cos(grd.u.latRad(:,2:end,:))...
                -u.val(:,1:end-1,:).*cos(grd.u.latRad(:,1:end-1,:)))...
                ./(grd.u.latRad(:,2:end,:)-grd.u.latRad(:,1:end-1,:));
            val0=(dv-du)./grd.p.latRad/earthRadius;
        else
            val0=0*val1;
        end

        val1=.25*val1(1:end-1,1:end-1,:)+.25*val1(2:end,1:end-1,:)...
            +.25*val1(2:end,1:end-1,:)+.25*val1(2:end,2:end,:); 
        val0=.25*val0(1:end-1,1:end-1,:)+.25*val0(2:end,1:end-1,:)...
            +.25*val0(2:end,1:end-1,:)+.25*val0(2:end,2:end,:);
        val=val1-val0;         
    else
        roms=get_roms_layer(roms_dir,fieldNames{iField},loc,'delta');
        if ~isempty(roms_ref)
            romsRef=get_roms_layer(roms_ref,fieldNames{iField},loc,'delta');
            val1=roms.val-romsRef.val;
        else
            val1=roms.val; 
        end
        val=val1(2:end-1,2:end-1,:);
    end

    
    %save to output
    out.(fieldNames{iField})=val; 
end

%% read 3D ROMS output


%3dfields
fieldNames={'temp','salt','u','v','vort'}; 

for iField=1:length(fieldNames)
    if ~ismember(fieldNames{iField},fieldOut); continue; end
    display(sprintf('reading %s',fieldNames{iField})); 
   
    
    %read
    loc=struct('t',t,'layer',layer);
    if strcmpi(fieldNames{iField},'u')
        roms=get_roms_layer(roms_dir,fieldNames{iField},loc,'delta');
        if ~isempty(roms_ref)
        romsRef=get_roms_layer(roms_ref,fieldNames{iField},loc,'delta');
        val1=roms.val-romsRef.val;
        else
            val1=roms.val;
        end
        val=.5*val1(1:end-1,2:end-1,:)+.5*val1(2:end,2:end-1,:);
    elseif strcmpi(fieldNames{iField},'v')
        roms=get_roms_layer(roms_dir,fieldNames{iField},loc,'delta');
        if ~isempty(roms_ref)
        romsRef=get_roms_layer(roms_ref,fieldNames{iField},loc,'delta');
        val1=roms.val-romsRef.val;
        else
            val1=roms.val; 
        end
        val=.5*val1(2:end-1,1:end-1,:)+.5*val1(2:end-1,2:end,:);
    elseif strcmpi(fieldNames{iField},'vort')
        u=get_roms_layer(roms_dir,'u',loc,'delta');
        v=get_roms_layer(roms_dir,'v',loc,'delta');
        
        dv=(v.val(2:end,:,:)-v.val(1:end-1,:,:))./...
            (grd.v.lonRad(2:end,:,:)-grd.v.lonRad(1:end-1,:,:)); 
        du=(u.val(:,2:end,:).*cos(grd.u.latRad(:,2:end,:))...
            -u.val(:,1:end-1,:).*cos(grd.u.latRad(:,1:end-1,:)))...
            ./(grd.u.latRad(:,2:end,:)-grd.u.latRad(:,1:end-1,:)); 
        val1=(dv-du)./grd.p.latRad/earthRadius; 
        
        if ~isempty(roms_ref)
            u=get_roms_layer(roms_ref,'u',loc,'delta');
            v=get_roms_layer(roms_ref,'v',loc,'delta');
            
            dv=(v.val(2:end,:,:)-v.val(1:end-1,:,:))./...
                (grd.v.lonRad(2:end,:,:)-grd.v.lonRad(1:end-1,:,:));
            du=(u.val(:,2:end,:).*cos(grd.u.latRad(:,2:end,:))...
                -u.val(:,1:end-1,:).*cos(grd.u.latRad(:,1:end-1,:)))...
                ./(grd.u.latRad(:,2:end,:)-grd.u.latRad(:,1:end-1,:));
            val0=(dv-du)./grd.p.latRad/earthRadius;
            
        else
            val0=0*val1;
        end
        
        val1=.25*val1(1:end-1,1:end-1,:)+.25*val1(2:end,1:end-1,:)...
            +.25*val1(2:end,1:end-1,:)+.25*val1(2:end,2:end,:); 
        val0=.25*val0(1:end-1,1:end-1,:)+.25*val0(2:end,1:end-1,:)...
            +.25*val0(2:end,1:end-1,:)+.25*val0(2:end,2:end,:);
        val=val1-val0; 
    else
        roms=get_roms_layer(roms_dir,fieldNames{iField},loc,'delta');
        if ~isempty(roms_ref)
        romsRef=get_roms_layer(roms_ref,fieldNames{iField},loc,'delta');
        val1=roms.val-romsRef.val;
        else
            val1=roms.val; 
        end
        val=val1(2:end-1,2:end-1,:);
    end


    %save to output
    out.(fieldNames{iField})=val; 
end

%% Add density

if ismember('rho',fieldOut)
    display('Calculating rho'); 
    p=sw_pres(squeeze(max(-grd.rho.z(2:end-1,2:end-1,layer),0)),squeeze(grd.rho.lat(2:end-1,2:end-1))); 
    p=repmat(p,[1,1,size(out.salt,3)]); 
    out.rho=sw_pden(out.salt,out.temp,p,0*p); 
    clear p; 
end

%% filter
if ~isempty(out_file{2})
    
    outf=out;
    
    if strcmpi(filter,'tidal')
        %u,v
        if isfield(out,'u') && isfield(out,'v')
            display('tidal filtering u,v');
            val=out.u+1i*out.v;
            dt=(t(2)-t(1))*24;
            for i2=1:size(val,2)
                for i1=1:size(val,1)
                    if all(isnan(squeeze(val(i1,i2,:)))); continue; end
                    tpar=t_tide(squeeze(val(i1,i2,:)),'interval',dt,...
                        'start time',t(1),'latitude',grd.rho.lat(i1+1,i2+1),'output','none');
                    tval=t_predic(t,tpar,'latitude',grd.rho.lat(i1+1,i2+1));
                    val(i1,i2,:)=val(i1,i2,:)-reshape(tval,1,1,[]);
                end
            end
            out.u=real(val); out.v=imag(val);
        end
        
    end
    
    %lanczo filter
    if strcmpi(filter,'lanczos')
        display('lanczos filter');
        
        fcut=2*fcut;
        for it=1:length(t)
            in=t>=t(it)-fwidth & t<=t(it)+fwidth;
            w=0.5*(1+cos(pi*(t(in)-t(it))/fwidth)).*sinc((t(in)-t(it))*fcut);
            w=w/sum(w);
            iB=find(in,1,'first'); iE=find(in,1,'last'); 
            for iField=1:length(fieldNames)
                if ~ismember(fieldNames{iField},fieldOut); continue; end
                outf.(fieldNames{iField})(:,:,it)=0; 
                for i0=iB:iE
                    outf.(fieldNames{iField})(:,:,it)=...
                        outf.(fieldNames{iField})(:,:,it)...
                        +w(i0-iB+1)*out.(fieldNames{iField})(:,:,i0); 
                end
            end
        end
    end
    
    %step
    if strcmpi(filter,'step')
        display('step filter');
        for it=1:length(t)
            in=t>=t(it)-fwidth & t<=t(it)+fwidth;
            w=ones(size(t(in))); 
            w=w/sum(w);
            iB=find(in,1,'first'); iE=find(in,1,'last');
            for iField=1:length(fieldNames)
                if ~ismember(fieldNames{iField},fieldOut); continue; end
                outf.(fieldNames{iField})(:,:,it)=...
                    w(1)*sum(out.(fieldNames{iField})(:,:,iB:iE),3);              
            end
        end
    end
    
end

    
%% save

out.t=t;
out.layer=layer;
out.source=roms_dir;
out.lon=grd.rho.lon(2:end-1,2:end-1);
out.lat=grd.rho.lat(2:end-1,2:end-1);
save(out_file{1},'-struct','out','-v7');
display(out_file{1}); 

if ~isempty(out_file{2})
    outf.t=t;
    outf.layer=layer;
    outf.source=roms_dir;
    outf.lon=grd.rho.lon(2:end-1,2:end-1);
    outf.lat=grd.rho.lat(2:end-1,2:end-1);
    save(out_file{2},'-struct','outf','-v7');
    display(out_file{2}); 
end


