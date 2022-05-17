%%%
%%% animVort.m
%%%
%%% Makes a movie of the vorticity.
%%%
%%% NOTE: Doesn't account for u/v gridpoint locations, and doesn't handle
%%% partial cells.
%%%

%%% Read experiment data
loadexp;
video_bool=0; % Make video or not. Set name on line 35

%%% Vertical grid spacing matrix
DZ = repmat(reshape(delR,[1 1 Nr]),[Nx Ny 1]);

%%% Diagnostic indix corresponding to instantaneous velocity
diagnum = length(diag_frequency);

%%% This needs to be set to ensure we are using the correct output
%%% frequency
diagfreq = diag_frequency(diagnum);

%%% Frequency of diagnostic output
dumpFreq = abs(diagfreq);
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters >= nIter0);
nDumps = length(dumpIters);

%%% Initialize movie
figure(8);
set(gcf,'Color','w');
M = moviein(nDumps);
if video_bool==1
    vid=VideoWriter('ref_tau1_Z20_vortEddies');
    vid.FrameRate=3;
    open(vid);
end

%%% Loop through iterations
for n=[1 11] %1:length(dumpIters)
% for n=240:nDumps
 
  tt(n) =  (dumpIters(n)-dumpIters(1))*deltaT/86400;
  tt(n);
  
  %%% Attempt to load either instantaneous velocities or their squares
  uvel = rdmdsWrapper(fullfile(exppath,'/results/UVEL_inst'),dumpIters(n)) ;      
  vvel = rdmdsWrapper(fullfile(exppath,'/results/VVEL_inst'),dumpIters(n)); 
  if (isempty(uvel) || isempty(vvel))   
    break;
  end
  
  %%% Plot the vorticity  
  [YY,XX] = meshgrid(yy,xx);  
  vort = zeros(Nx,Ny);
  zlev = 8;
  vort(:,2:Ny) = - (uvel(:,2:Ny,zlev)-uvel(:,1:Ny-1,zlev))/delY(1);
  vort = vort + (vvel([2:Nx 1],:,zlev)-vvel(:,:,zlev))/delX(1);  
  ff = f0*ones(Nx,Ny);
  vort=vort./abs(ff); 
  pcolor(XX/1000,YY/1000,vort)
  % contourf(XX/1000,YY/1000,vort./abs(ff),[0.2 0.5 1]);
  shading interp;
  %   contourf(XX/1000,YY/1000,vort./abs(ff),-.4:0.005:.4,'EdgeColor','None');
  %   hold on;
%   contour(XX/1000,YY/1000,-bathy,500:500:3500,'EdgeColor',[0.5 0.5 0.5]);
%   hold off;
  colorbar;
  colormap redblue;
  caxis([-1 1]);
  set(gca,'FontSize',12);
  xlabel('x (km)');
  ylabel('y (km)');
%   set(gca,'XTick',[-800:400:800]);
%   set(gca,'YTick',[0:400:2000]);
  title(['\zeta/|f| at t= ',num2str(round(tt(n)),'%3d'),' days']);
  M(n) = getframe(gcf);
  if video_bool==1
      writeVideo(vid,M(n));
  end
end
if video_bool==1
    close(vid);
end
