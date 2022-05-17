%%%
%%% anim.m
%%%
%%% Reads diagnostic output from MITgcm and makes a movie of it. Requires
%%% that the experiment be generated using the 'newexp' package.
%%%

%%% NOTE: Doesn't account for u/v gridpoint locations

%%% Set true if plotting on a Mac
mac_plots = 0;

%%% Read experiment data
loadexp;

%%% Select diagnostic variable to animate
diagnum = 6;
outfname = diag_fileNames{1,diagnum};

%%% Data index in the output data files
outfidx = 1;

%%% If set true, plots a top-down view of the field in a given layer.
%%% Otherwise plots a side-on view of the zonally-averaged field.
xyplot = 1;

%%% Vertical layer index to use for top-down plots
xylayer = 20;

%%% Set true to plot the field in the lowest active cell at each horizontal
%%% location
botplot = 0;

%%% Set true for a zonal average
yzavg = 0;

%%% Layer to plot in the y/z plane
yzlayer = 1;

%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics.
dumpFreq = abs(diag_frequency(diagnum));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);

%%% Mesh grids for plotting
hFac = hFacC;
kmax = zeros(Nx,Ny);
if (xyplot)
  [YY,XX] = meshgrid(yy/1000,xx/1000);
  kmax = sum(ceil(hFac),3);
  kmax(kmax==0) = 1;
else  
  %%% Create mesh grid with vertical positions adjusted to sit on the bottom
  %%% topography and at the surface
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
handle = figure(9);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
set(gcf,'color','w');
M = moviein(nDumps);

Amean = [];
tdays = [];
for n=1:length(dumpIters)
% for n=50:length(dumpIters)
  
  t = dumpIters(n)*deltaT/86400;
  
  if (n > 1)
    Aprev = A;
  end
%   A = rdmdsWrapper(fullfile(exppath,'results',outfname),0);          
  A = rdmdsWrapper(fullfile(exppath,'results',outfname),dumpIters(n));          
  if (isempty(A))
    error(['Ran out of data at t=,',num2str(tdays(n)),' years']);
  end      
  
  max(max(max(abs(A))))
 
  tdays(n) = t;
  
  %%% x/y plot
  if (xyplot)
    
    if (botplot)      
      FF = zeros(Nx,Ny);
      for i=1:Nx
        for j=1:Ny
          FF(i,j) = A(i,j,kmax(i,j));
%           FF(i,j) = A(i,j,kmax(i,j))*hFac(i,j,kmax(i,j));
        end
      end
    else
      FF = squeeze(A(:,:,xylayer,outfidx));        
    end
    
    FF(hFacC(:,:,xylayer)==0) = NaN;
%     contourf(XX,YY,FF,100,'EdgeColor','None');  
    pcolor(XX,YY,FF);
    shading interp;
    xlabel('x (km)');
    ylabel('y (km)');
    colormap jet;
%     caxis([-0.4 0])
%     caxis([-6 5]);
    
  %%% y/z zonally-averaged plot
  else
    
    if (yzavg)
      Ayz = squeeze(nanmean(A(:,:,:)));    
    else
      Ayz = squeeze(A(yzlayer,:,:,outfidx));
    end
    
    jrange = 1:Ny;
%     [C h] = contourf(YY(jrange,:),ZZ(jrange,:)/1000,Ayz(jrange,:),200,'EdgeColor','None');              
    pcolor(YY(jrange,:),ZZ(jrange,:)/1000,Ayz(jrange,:));
    shading interp;
    hold on;
%     [C,h]=contour(YY(jrange,:),ZZ(jrange,:)/1000,Ayz(jrange,:),10,'EdgeColor','k');
%     [C,h]=contour(YY(jrange,:),ZZ(jrange,:)/1000,Ayz(jrange,:),[-2:0.5:12],'EdgeColor','k');
%     [C,h]=contour(YY(jrange,:),ZZ(jrange,:)/1000,Ayz(jrange,:),[-0.1:0.005:0.1],'EdgeColor','k');    
    hold off;
    colormap jet;
    xlabel('Offshore $y$ (km)','interpreter','latex');
    ylabel('Height $z$ (km)','interpreter','latex');  
  
  end
    
  %%% Finish the plot
  handle=colorbar;
  set(handle,'FontSize',fontsize);
  set(gca,'FontSize',fontsize);
  title(['$t=',num2str(tdays(n),'%.1f'),'$ days'],'interpreter','latex');  
  set(gca,'Position',plotloc);
%   annotation('textbox',[0.85 0.05 0.25 0.05],'String','$\overline{\theta}$ ($^\circ$C)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
  M(n) = getframe(gcf);
  
end