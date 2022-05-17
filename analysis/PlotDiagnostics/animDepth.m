%%%
%%% anim.m
%%%
%%% Reads diagnostic output from MITgcm and makes a movie of it. Requires
%%% that the experiment be generated using the 'newexp' package.
%%%

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
xylayer = 1;


%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics.
dumpFreq = abs(diag_frequency(diagnum));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);

%%% Mesh grids for plotting
hFac = hFacC;
kmax = zeros(Nx,Ny);
[YY,XX] = meshgrid(yy/1000,xx/1000);
kmax = sum(ceil(hFac),3);
kmax(kmax==0) = 1;


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
if make_movie == 1
    vid=VideoWriter('Ref_Z'); 
    vid.FrameRate=3;
    open(vid);
end
Amean = [];
tdays = [];
for n=[1 3 10 28] %1:length(dumpIters)-2
  
  t = dumpIters(n)*deltaT/86400;
  
  if (n > 1)
    Aprev = A;
  end
  
  A = rdmdsWrapper(fullfile(exppath,'results',outfname),dumpIters(n));          
  if (isempty(A))
    error(['Ran out of data at t=,',num2str(tdays(n)),' years']);
  end      
  
  if (n > 1)
    ['diff = ',num2str(max(max(max(abs(A-Aprev)))))];
  end     
 
  tdays(n) = t;
  
  %%% x/y plot
  for xylayer=1:length(zz)/2
      FF = squeeze(A(:,:,xylayer,outfidx));
      FF(hFacC(:,:,xylayer)==0) = NaN;
      % contourf(XX,YY,FF,100,'EdgeColor','None');
      pcolor(XX,YY,FF);
      shading interp;
      xlabel('x (km)');
      ylabel('y (km)');
      colormap jet;
      caxis([34.2 34.6])
      
      %%% Finish the plot
      handle=colorbar;
      set(handle,'FontSize',fontsize);
      set(gca,'FontSize',fontsize);
      % title(['$t=',num2str(tdays(n),'%.1f'),'$ days'],'interpreter','latex');
      title(['$z=',num2str(zz(xylayer),'%.1f'),'$ days'],'interpreter','latex');
      set(gca,'Position',plotloc);
      %   annotation('textbox',[0.85 0.05 0.25 0.05],'String','$\overline{\theta}$ ($^\circ$C)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
      M(n) = getframe(gcf);
      if make_movie == 1
          writeVideo(vid,M(n))
      end
  end
end

if make_movie == 1
    close(vid)
end

save(expname)