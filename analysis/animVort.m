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

%%% Vertical level to plot
zlev = 20;

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

%%% Loop through iterations
for n=1:nDumps
% for n=240:nDumps
 
  tt(n) =  (dumpIters(n)-dumpIters(1))*deltaT/86400;
  tt(n)
  
  %%% Attempt to load either instantaneous velocities or their squares
  uvel = rdmdsWrapper(fullfile(exppath,'/results/UVEL_inst'),dumpIters(n));      
  vvel = rdmdsWrapper(fullfile(exppath,'/results/VVEL_inst'),dumpIters(n)); 
  if (isempty(uvel) || isempty(vvel))   
    break;
  end
  
  %%% Plot the vorticity  
  [YY,XX] = meshgrid(yy,xx);  
  vort = zeros(Nx,Ny);  
  vort(:,2:Ny) = - (uvel(:,2:Ny,zlev)-uvel(:,1:Ny-1,zlev))/delY(1);
  vort = vort + (vvel([2:Nx 1],:,zlev)-vvel(:,:,zlev))/delX(1);
%   ubt = sum(uvel.*DZ.*hFacW,3) ./ sum(DZ.*hFacW,3);
%   vbt = sum(vvel.*DZ.*hFacS,3) ./ sum(DZ.*hFacS,3);
%   vort(:,2:Ny) = - (ubt(:,2:Ny)-ubt(:,1:Ny-1))/delY(1);
%   vort = vort + (vbt([2:Nx 1],:)-vbt(:,:))/delX(1);  
  ff = f0*ones(Nx,Ny);
  pcolor(XX/1000,YY/1000,vort./abs(ff));
  shading interp;
%   contourf(XX/1000,YY/1000,vort./abs(ff),-.4:0.005:.4,'EdgeColor','None');
%   hold on;
%   contour(XX/1000,YY/1000,-bathy,500:500:3500,'EdgeColor',[0.5 0.5 0.5]);
%   hold off;
  colorbar;
  colormap redblue;
%   caxis([-0.4 0.4]);
%   caxis([-0.1 0.1]);
  caxis([-3 3]);
  set(gca,'FontSize',12);
  xlabel('x (km)');
  ylabel('y (km)');
%   set(gca,'XTick',[-800:400:800]);
%   set(gca,'YTick',[0:400:2000]);
  title(['\zeta/|f| at t= ',num2str(round(tt(n)),'%3d'),' days']);
  M(n) = getframe(gcf);  
  
end