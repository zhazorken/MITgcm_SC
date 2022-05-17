%%%
%%% animV2.m
%%%
%%% Reads diagnostic output from MITgcm and makes a movie of it. Requires
%%% that the experiment be generated using the 'newexp' package.
%%%
%%% Uses subplots to animate multiple runs (rows) on the same figure for:
%%% Column 1 : density along the lead, meridionally averaged
%%% Column 2 : vorticity at mixed layer depth calculated by the function
%%% HmixMIT
%%% Column 3 : vorticity at the surface

%%% NOTE: Doesn't account for u/v gridpoint locations

%%% Add paths
addpath /data3/kayliec23/MITgcm_LC/GSW;
addpath /data3/kayliec23/MITgcm_LC/GSW/html;
addpath /data3/kayliec23/MITgcm_LC/GSW/library;
addpath /data3/kayliec23/MITgcm_LC/GSW/pdf;

%%% Set true to make movie
make_movie = 0; % Set name on line 119

%%% Set true if plotting on a Mac
mac_plots = 0;

%%% Choose files to animate
outlist{1}='/data3/kayliec23/MITgcm_LC/experiments/Ref_Hpyc';
expName{1}='Hpyc50';
outlist{2}='/data3/kayliec23/MITgcm_LC/experiments/Ref_Hpyc';
expName{2}='Hpyc400';
outlist{3}='/data3/kayliec23/MITgcm_LC/experiments/Ref_Tair';
expName{3}='Tair10';
outlist{4}='/data3/kayliec23/MITgcm_LC/experiments/Ref_Tair';
expName{4}='Tair40';
outlist{5}='/data3/kayliec23/MITgcm_LC/experiments/Ref_H0';
expName{5}='H010';
outlist{6}='/data3/kayliec23/MITgcm_LC/experiments/Ref_H0';
expName{6}='H0200';
outlist{7}='/data3/kayliec23/MITgcm_LC/experiments/Ref_W';
expName{7}='W1';
outlist{8}='/data3/kayliec23/MITgcm_LC/experiments/Ref_W';
expName{8}='W5';
outlist{9}='/data3/kayliec23/MITgcm_LC/experiments/Ref_Tmax';
expName{9}='Tmax0';
outlist{10}='/data3/kayliec23/MITgcm_LC/experiments/Ref_Tmax';
expName{10}='Tmax1';
outlist{11}='/data3/kayliec23/MITgcm_LC/experiments/Ref_drag';
expName{11}='v001';
outlist{12}='/data3/kayliec23/MITgcm_LC/experiments/Ref_drag';
expName{12}='v01_redo';
outlist{13}='/data3/kayliec23/MITgcm_LC/experiments/Ref_dS';
expName{13}='dS_min';
outlist{14}='/data3/kayliec23/MITgcm_LC/experiments/Ref_dS';
expName{14}='dS_max';

acrossLead_bool=1;
alongLead_bool=0;
vortSurf_bool=0;
vortML_bool=0;
dens_surf=0;
densML=0;

subplotRows=2; %length(outlist);
subplotCols=5; %length(outlist);
numDays=30; % How many days to simulate
Hpyc=[50 400];
count=0;
botplot=0;
wstart=189;
wend=211;

% Setting up movie
if make_movie == 1
    %figure(9)
    vid=VideoWriter('Ref_Hpyc');
    vid.FrameRate=3;
    open(vid);
end

for n=[1 numDays]
        plotPos=0;
        count=count+1;
    %%% Loop over files
    for f=6 %1:length(outlist)
        basedir=outlist{f};
        expdir=basedir;
        expname=expName{f};
        loadexp;
        
        
        %%% Select diagnostic variable to animate
        diagnum = 6;
        density_bool=1; % Set true to plot density
        outfname = diag_fileNames{1,diagnum};
        %%% Data index in the output data files
        outfidx = 1;
        
        
        %%% If set true, plots a top-down view of the field in a given layer.
        %%% Otherwise plots a side-on view of the zonally-averaged field.
        xyplot = 0;
        yzplot = 1; % plot across the lead
        xzplot = 0; % plot along the lead
        
        %%% Vertical layer index to use for top-down plots
        yzlayer = length(yy)/2; % across
        xzlayer = length(xx)/2; % along
        
        %%% Set true to plot the field in the lowest active cell at each horizontal
        %%% location
        botplot = 0;
        
        %%% Set true for a meridional/zonal average
        xzavg = 1; % along
        yzavg = 1; % across
        
        %%% Frequency of diagnostic output - should match that specified in
        %%% data.diagnostics.
        dumpFreq = abs(diag_frequency(diagnum));
        nDumps = round(nTimeSteps*deltaT/dumpFreq);
        dumpIters = round((1:nDumps)*dumpFreq/deltaT);
        dumpIters = dumpIters(dumpIters > nIter0);
        
        
        %%% Mesh grids for plotting
        hFac = hFacC;
        kmax = zeros(Nx,Ny);
        if xyplot==1
            [YY,XX] = meshgrid(yy/1000,xx/1000);
            kmax = sum(ceil(hFac),3);
            kmax(kmax==0) = 1;
        end
            
        if xzplot==1
            [ZZ,XX] = meshgrid(zz,xx/1000);
            for i=1:Nx
                if (xzavg)
                    hFacC_row = squeeze(hFacC(i,:,:));
                    hFacC_row = max(hFacC_row,[],1);
                else
                    hFacC_row = squeeze(hFacC(i,xzlayer,:))';
                end
                kmax = length(hFacC_row(hFacC_row>0));
                zz_botface = -sum(hFacC_row.*delR);
                ZZ(i,1) = 0;
                if (kmax>0)
                    ZZ(i,kmax) = zz_botface;
                end
            end
        end
        if yzplot==1
            [ZZ,YY] = meshgrid(zz,yy/1000);
            for j=1:Ny
                if (yzavg)
                    hFacC_col = squeeze(hFacC(:,j,:));
                    hFacC_col = max(hFacC_col,[],1);
                else
                    hFacC_col = squeeze(hFacC(yzlayer,j,:))';
                end
                kmax = length(hFacC_col(hFacC_col>0));
                zz_botface = -sum(hFacC_col.*delR);
                ZZ(j,1) = 0;
                if (kmax>0)
                    ZZ(j,kmax) = zz_botface;
                end
            end
        end
        
        %%% Plotting options
        scrsz = get(0,'ScreenSize');
        fontsize = 26;
        if (mac_plots)
            framepos = [0 scrsz(4)/2 scrsz(3)/1.3 scrsz(4)];
            plotloc = [0.17 0.3 0.62 0.7];
        else
            plotloc = [0.11 0.15 0.76 0.73];
            framepos = [327    80   941   885];
        end
        
        %%% Set up the figure
        if f==1 && n==1
            handle = figure(9);
            set(handle,'Position',framepos);
            clf;
            axes('FontSize',fontsize);
            set(gcf,'color','w');
        end
        % M = moviein(nDumps);
        Amean = [];
        tdays = [];
        
        
        t = dumpIters(n)*deltaT/86400;
        
% % %         if (n > 1)
% % %             Aprev = A;
% % %         end
% % %         
% % %         A = rdmdsWrapper(fullfile(exppath,'results',outfname),dumpIters(n));
% % %         if (isempty(A))
% % %             error(['Ran out of data at t=,',num2str(tdays(n)),' years']);
% % %         end
% % %         
% % %         if (n > 1)
% % %             ['diff = ',num2str(max(max(max(abs(A-Aprev)))))];
% % %         end
        
        % Across-Lead Velocity
        outfname_velU = diag_fileNames{1,2};
        if (n > 1)
            velUAprev = velU;
        end
        
        velU = rdmdsWrapper(fullfile(exppath,'results',outfname_velU),dumpIters(n));
        if (isempty(velU))
            error(['Ran out of data at t=,',num2str(tdays(n)),' years']);
        end
        
        if (n > 1)
            ['diff = ',num2str(max(max(max(abs(velU-velUAprev)))))];
        end
        
        % Salinity
        outfname_sal = diag_fileNames{1,6};
        if (n > 1)
            salAprev = sal;
        end
        
        sal = rdmdsWrapper(fullfile(exppath,'results',outfname_sal),dumpIters(n));
        if (isempty(sal))
            error(['Ran out of data at t=,',num2str(tdays(n)),' years']);
        end
        
        if (n > 1)
            ['diff = ',num2str(max(max(max(abs(sal-salAprev)))))];
        end
        
        % Temperature
        outfname_temp = diag_fileNames{1,5};
        if (n > 1)
            tempAprev = temp;
        end
        
        temp = rdmdsWrapper(fullfile(exppath,'results',outfname_temp),dumpIters(n));
        if (isempty(temp))
            error(['Ran out of data at t=,',num2str(tdays(n)),' years']);
        end
        
        if (n > 1)
            ['diff = ',num2str(max(max(max(abs(temp-tempAprev)))))];
        end
        
        % Density
        
        pdens=gsw_rho(sal,temp,0);
        if density_bool==1
            A=pdens;
        end
        
        press_=zeros(1, 1,size(sal,3));
        press_(1,1,:)=linspace(0,495,100);
        press=repmat(press_,[400,400,1]);
        
        dens=gsw_rho(sal,temp,press);        
        
        tdays(n) = t;
       
        
        %%%%%%%%%%%% Plotting meridional and zonal plots %%%%%%%%%%%%%
        
        %%% x/z meriodionally-averaged plot     
        
        %     [C h] = contourf(YY(jrange,:),ZZ(jrange,:)/1000,Ayz(jrange,:),200,'EdgeColor','None');
        % if n==4 || n==30
        
        if alongLead_bool==1
            Axz = squeeze(nanmean(A(:,wstart:wend,:),2));
            jrange = 1:Nx;
            [ZZ XX]=meshgrid(zz,xx);
            plotPos=plotPos+1;
            sp1=subplot(subplotRows,subplotCols,plotPos);
            pcolor(XX(jrange,:),ZZ(jrange,:)/1000,Axz(jrange,:));
            shading interp;
            colormap(sp1,pmkmp(256,'IsoAZ180'))
            colorbar
            xlabel('$y$ (km)','interpreter','latex');
            ylabel('$z$ (km)','interpreter','latex');
            hold on;
        end
        %%% y/z zonally-averaged plot

        % [C h] = contourf(YY(jrange,:),ZZ(jrange,:)/1000,Ayz(jrange,:),200,'EdgeColor','None');
        %if n==4 || n==30
        if acrossLead_bool==1
            figure(2)
            [ZZ YY]=meshgrid(zz,yy);
            Ayz = squeeze(nanmean(A(:,:,:)));
            jrange = 1:Ny;
            plotPos=plotPos+1;
            %sp2=subplot(subplotRows,subplotCols,plotPos);
            pcolor(YY(jrange,:)/1000,ZZ(jrange,:),Ayz(jrange,:)); 
            colormap(pmkmp(256,'Swtth'));
            %colormap(sp2,pmkmp(256,'Swtth'))
            shading interp
            ylim([-350 0])
            if f==length(outlist)
            xlabel('$y$ (km)','interpreter','latex');
            end
            if f==1
            ylabel('$z$ (m)','interpreter','latex');
            colorbar
            end
            hold on;
            %plot(Hm_line, 'LineWidth',1.5,'Color','k')
            %pbaspect([2.5 1 1])
            %caxis([1027.42 1027.6])
            colorbar
        end
        
        %%%%%%%%%%%% Plotting top-down plots %%%%%%%%%%%%%
        
        %%% Density %%%
        
        
        [YY XX ZZ]=meshgrid(yy,xx,zz);
        % plot density at surface
        if dens_surf==1
            zlev=2;
            if xyplot==1
                if (botplot)
                    FF = zeros(Nx,Ny);
                    for i=1:Nx
                        for j=1:Ny
                            FF(i,j) = pdens(i,j,kmax(i,j));
                        end
                    end
                else
                    FF = squeeze(pdens(:,:,zlev));
                end
            end
            FF(hFacC(:,:,zlev)==0) = NaN;
            plotPos=plotPos+1;
            sp5=subplot(subplotRows,subplotCols,plotPos);
            pcolor(XX,YY,FF);
            shading interp;
            colormap(sp5,pmkmp(256,'Swtth'))
            caxis([1027.42 1027.6])
            colorbar
       end
            
            % plot density at Hm
            if densML==1
                zlev=Hm_lev(n);
                if xyplot==1
                    if (botplot)
                        FF = zeros(Nx,Ny);
                        for i=1:Nx
                            for j=1:Ny
                                FF(i,j) = pdens(i,j,kmax(i,j));
                            end
                        end
                    else
                        FF = squeeze(pdens(:,:,zlev));
                    end
                end
                
                FF(hFacC(:,:,zlev)==0) = NaN;
                isovalue=pdens(Nx/2,Ny/2,Hm_lev(n)); %1027.48; %nanmean(nanmean(pdens(:,wstart:wend,zlev)));
                %plotPos=plotPos+1;
                %sp5=subplot(subplotRows,subplotCols,plotPos);
                %pcolor(XX,YY,FF);
                fv=isosurface(XX/1000,YY/1000,ZZ,pdens,isovalue);
                p=patch(fv);
                p.FaceColor=[79 66 181]/225;
                p.EdgeColor='none';
                %shading interp;
                %colormap(sp5,pmkmp(256,'Swtth'))
                %caxis([1027.46 1027.55])
                %colorbar
                hold on;
            end
            
        
        %%% Vorticity %%%
               
        %%% Vertical grid spacing matrix
        DZ = repmat(reshape(delR,[1 1 Nr]),[Nx Ny 1]);
        
        tt(n) =  (dumpIters(n)-dumpIters(1))*deltaT/86400;
        tt(n);
        
        %%% Attempt to load either instantaneous velocities or their squares
        uvel = rdmdsWrapper(fullfile(exppath,'/results/UVEL_inst'),dumpIters(n)) ;
        vvel = rdmdsWrapper(fullfile(exppath,'/results/VVEL_inst'),dumpIters(n));
        
        if count==1
            v = zeros(length(outlist),numDays,length(yy),length(zz));
        end
        v(f,n,:,:) = mean(vvel);
        
        if (isempty(uvel) || isempty(vvel))
            break;
        end
        
        %%% Plot the vorticity at Hmix
        
        [YY,XX] = meshgrid(yy,xx);
        vort = zeros(Nx,Ny);
        
%         vort(:,2:Ny) = - (uvel(:,2:Ny,zlev)-uvel(:,1:Ny-1,zlev))/delY(1);
%         vort = vort + (vvel([2:Nx 1],:,zlev)-vvel(:,:,zlev))/delX(1);
        ff = f0*ones(Nx,Ny);
        
        
        % Plot vorticity at surface
        if vortSurf_bool==1
            zlev=2;
            vort(:,2:Ny) = - (uvel(:,2:Ny,zlev)-uvel(:,1:Ny-1,zlev))/delY(1);
            vort = vort + (vvel([2:Nx 1],:,zlev)-vvel(:,:,zlev))/delX(1);
            vort=vort./abs(ff);
            %plotPos=plotPos+1;
            %sp3=subplot(subplotRows,subplotCols,plotPos);
            %pcolor(XX/1000,YY/1000,vort./abs(ff))
            p=surf(XX/1000,YY/1000,0*XX,vort);
            p.FaceColor='texturemap';
            colormap redblue
            shading interp
            % colormap (sp3,redblue)
            colorbar
            p.EdgeColor='none';
            caxis([-1 1]);
            set(gca,'FontSize',12);
            %if n==30
                xlabel('x (km)');
            %end
            ylabel('y (km)');
            % title(['\zeta/|f| at t= ',num2str(round(tt(n)),'%3d'),' days']);
            hold on
            alpha(p,0.5);
        end
       
        
        if vortML_bool==1
            %figure(6)
            zlev=Hm_lev(n);
            vort(:,2:Ny) = - (uvel(:,2:Ny,zlev)-uvel(:,1:Ny-1,zlev))/delY(1);
            vort = vort + (vvel([2:Nx 1],:,zlev)-vvel(:,:,zlev))/delX(1);
            vort=vort./abs(ff);
            plotPos=plotPos+1;
            sp3=subplot(subplotRows,subplotCols,plotPos);
            pcolor(XX/1000,YY/1000,vort)
            shading interp
            colormap(sp3,redblue)
            caxis([-1 1]);
            set(gca,'FontSize',12);
            if f==length(outlist)
                xlabel('x (km)');
            end
            if f==1
            ylabel('y (km)');
            title(['\zeta/|f| at t= ',num2str(round(tt(n)),'%3d'),' days']);
            colorbar
            end
            hold on
            pbaspect([1 1 1])
        end
        
%         
% %         % Eddy Diffusivity
%         dens = pdens; 
%         dens_v = dens -repmat(mean(dens,1),[Nx 1 1]);
%         vvel = vvel - repmat(mean(vvel,1),[Nx 1 1]);
%         dens_v = 0.5 * (dens_v(:,1:Ny,:)+dens_v(:,[Ny 1:Ny-1],:));
%         num=mean(mean(vvel.*dens_v,3),1);    
%         denom=mean(mean((dens(:,1:Ny,:)-dens(:,[Ny 1:Ny-1],:))/delY(1),3),1);
%         kappa=-num./denom ;
%         
%         % Saving data
%         if kappaFilter==1
%             for i=wstart-49:wend+49
%                 kappa(i)=nan;
%             end
%         end
%         
%         results(f).data(count).eddyDiff=kappa;
%            
%         plotPos=plotPos+1;
%         sp4=subplot(subplotRows,subplotCols,plotPos);
%         
        
        
        
        
        
            
            %         %%% Finish the plot
            %         handle=colorbar;
            %         set(handle,'FontSize',fontsize);
            %         set(gca,'FontSize',fontsize);
            %         title(['$t=',num2str(tdays(n),'%.1f'),'$ days'],'interpreter','latex');
            %         set(gca,'Position',plotloc);
            %         %   annotation('textbox',[0.85 0.05 0.25 0.05],'String','$\overline{\theta}$ ($^\circ$C)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
            %    figlabel = strcat(['t= ',num2str(round(tt(n)),'%3d'),' days']);
            %    delete(findall(gcf,'tag','annot'));
%             if n==2
%                 pos=[0.01 0.82 0.06 0.01];
%             elseif n==4
%                 if f==1
%                     pos=[0.01 0.95 0.1 0.1];
%                 elseif f==2
%                     pos=[0.01 0.7 0.1 0.1];
%                 end
%                 %figlabel = strcat([{'H_{pyc}=',num2str(Hpyc(f)),' m'},{'t= ',num2str(round(tt(n)),'%3d'),' days'}]);
%                 figlabel = strcat(['t= ',num2str(round(tt(n)),'%3d'),' days']);
%             elseif n==6
%                 pos=[0.01 0.58 0.1 0.01];
%             elseif n==11
%                 pos=[0.01 0.43 0.1 0.01];
%             elseif n==31
%                 if f==1
%                     pos=[0.01 0.45 0.1 0.1];
%                 elseif f==2
%                     pos=[0.01 0.2 0.1 0.1];
%                 end
%                 figlabel = strcat(['t= ',num2str(round(tt(n)),'%3d'),' days']);
%             end
            %    if n>1
            %        annotation('textbox',pos,'String',figlabel,'Fontsize',12,'LineStyle','none','Tag','annot')
            %    end
            %figlabel = strcat(['t= ',num2str(round(tt(n)),'%3d'),' days']);
            %delete(findall(gcf,'tag','annot'));
            %annotation('textbox',[0.5 0.99 0.1 0.01],'String',figlabel,'Fontsize',16,'LineStyle','none','Tag','annot')
            %M(n) = getframe(gcf);
            %if make_movie == 1
            %   writeVideo(vid,M(n))
    end
end

if make_movie == 1
    close(vid)
end
save(expname)
