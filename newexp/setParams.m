%%%
%%% setParams.m
%%%
%%% Sets basic MITgcm parameters plus parameters for included packages, and
%%% writes out the appropriate input files.,
%%%
function nTimeSteps = setParams (inputpath,codepath,listterm,Nx,Ny,Nr)  

  %%% Load EOS utilities
  addpath ~/UCLA/Utilities/GSW;
  addpath ~/UCLA/Utilities/GSW/html;
  addpath ~/UCLA/Utilities/GSW/library;
  addpath ~/UCLA/Utilities/GSW/pdf;
    
  %%% TODO ADD SEA ICE AND ATMOSPHERIC FORCING; REMOVE SW, LW, SENS HEAT FLUXES      
  %%% TODO add sponges for waves?
  
  %%%%%%%%%%%%%%%%%%
  %%%%% SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%      
  
  %%% If set true, plots of prescribed variables will be shown
  showplots = true;      
  fignum = 1;
  
  %%% Data format parameters
  ieee='b';
  prec='real*8';
  realdigits = 8;
  realfmt=['%.',num2str(realdigits),'e'];
  
  %%% Get parameter type definitions
  paramTypes;
  
  %%% To store parameter names and values
  parm01 = parmlist;
  parm02 = parmlist;
  parm03 = parmlist;
  parm04 = parmlist;
  parm05 = parmlist;
  PARM={parm01,parm02,parm03,parm04,parm05}; 
  
  %%% Seconds in one hour
  t1min = 60;
  %%% Seconds in one hour
  t1hour = 60*t1min;
  %%% hours in one day
  t1day = 24*t1hour;
  %%% Seconds in 1 year
  t1year = 365*t1day;  
  %%% Metres in one kilometre
  m1km = 1000; 
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% SIMULATION CONTROL PARAMETERS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  nonHydrostatic = true; %%% Whether to run in nonhydrostatic mode
  use_seaIce = false; %%% Whether to run with sea ice (not yet implemented)
  use_3D = false; %%% Whether to run a 3D vs 2D simulation
  Ypoly = 0; %%% Latitudinal location of polynya
  Wpoly = 5*m1km; %%% Latitudinal width of polynya
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% FIXED PARAMETER VALUES %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  
  simTime = 30*t1day; %%% Simulation time  
  nIter0 = 0; %%% Initial iteration 
  Lx = 50*m1km; %%% Domain size in x 
  Ly = 50*m1km; %%% Domain size in y   
  H = 500; %%% Domain size in z 
  g = 9.81; %%% Gravity
  Omega = 2*pi*366/365/86400;  
  lat0 = -70; %%% Latitude at southern boundary
  f0 = 2*Omega*sind(lat0); %%% Coriolis parameter      
  rho0 = 1000;  
  
  %%% Diffusion parameters
  viscAh = 0; %%% Horizontal viscosity    
  viscA4 = 0; %%% Biharmonic viscosity
  viscAhGrid = 0; %%% Grid-dependent viscosity
%   viscA4Grid = 0.2; %%% Grid-dependent biharmonic viscosity    
%   viscC4smag = 0; %%% Smagorinsky biharmonic viscosity
%   diffK4Tgrid = 0.2; %%% Grid-dependent biharmonic diffusivity
  viscA4Grid = 0; %%% Grid-dependent biharmonic viscosity    
  viscC4smag = 4; %%% Smagorinsky biharmonic viscosity
  diffK4Tgrid = 0.1; %%% Grid-dependent biharmonic diffusivity
  viscAr = 1e-4; %%% Vertical viscosity
  diffKhT = 0; %%% Horizontal temp diffusion
  diffKrT = 0; %%% Vertical temp diffusion     
  
  %%% Parameters related to periodic forcing
  periodicExternalForcing = true;
  externForcingCycle = 2*simTime;  
  if (~periodicExternalForcing)
    externForcingPeriod = externForcingCycle;
    nForcingPeriods = 1;
  else
    externForcingPeriod = 0.01*t1day;
    nForcingPeriods = externForcingCycle/externForcingPeriod;
  end
  
  
  %%% PARM01
  %%% momentum scheme
  parm01.addParm('vectorInvariantMomentum',true,PARM_BOOL);
  %%% viscosity  
  parm01.addParm('viscAr',viscAr,PARM_REAL);
  parm01.addParm('viscA4',viscA4,PARM_REAL);
  parm01.addParm('viscAh',viscAh,PARM_REAL);
  parm01.addParm('viscA4Grid',viscA4Grid,PARM_REAL);
  parm01.addParm('viscAhGrid',viscAhGrid,PARM_REAL);
  parm01.addParm('viscA4GridMax',0.5,PARM_REAL);
  parm01.addParm('viscAhGridMax',1,PARM_REAL);
  parm01.addParm('useAreaViscLength',false,PARM_BOOL);
  parm01.addParm('useFullLeith',true,PARM_BOOL);
  parm01.addParm('viscC4leith',0,PARM_REAL);
  parm01.addParm('viscC4leithD',0,PARM_REAL);      
  parm01.addParm('viscC2leith',0,PARM_REAL);
  parm01.addParm('viscC2leithD',0,PARM_REAL); 
  parm01.addParm('viscC2smag',0,PARM_REAL); 
  parm01.addParm('viscC4smag',viscC4smag,PARM_REAL); 
  
  %%% diffusivity
  parm01.addParm('diffKrT',diffKrT,PARM_REAL);
  parm01.addParm('diffKhT',diffKhT,PARM_REAL);  
  %%% advection schemes
  parm01.addParm('tempAdvScheme',33,PARM_INT);
  parm01.addParm('saltAdvScheme',33,PARM_INT);
%   parm01.addParm('tempAdvScheme',7,PARM_INT);
%   parm01.addParm('saltAdvScheme',7,PARM_INT);
  parm01.addParm('multiDimAdvection',true,PARM_BOOL);
  parm01.addParm('tempStepping',true,PARM_BOOL);
  parm01.addParm('saltStepping',true,PARM_BOOL);
  parm01.addParm('staggerTimeStep',true,PARM_BOOL);
%   parm01.addParm('convertFW2Salt',-1,PARM_REAL);
  %%% equation of state
  parm01.addParm('eosType','JMD95Z',PARM_STR);   
  %%% boundary conditions
  parm01.addParm('no_slip_sides',false,PARM_BOOL);
  parm01.addParm('no_slip_bottom',false,PARM_BOOL);
  parm01.addParm('bottomDragLinear',0,PARM_REAL);
  parm01.addParm('bottomDragQuadratic',0,PARM_REAL);
  %%% physical parameters
  parm01.addParm('f0',f0,PARM_REAL);
  parm01.addParm('gravity',g,PARM_REAL);
  parm01.addParm('rhonil',rho0,PARM_REAL);
  parm01.addParm('rhoConst',rho0,PARM_REAL);
 
  %%% implicit diffusion and convective adjustment  
  parm01.addParm('ivdc_kappa',0,PARM_REAL);
  parm01.addParm('implicitDiffusion',~nonHydrostatic,PARM_BOOL);
  parm01.addParm('implicitViscosity',~nonHydrostatic,PARM_BOOL);
  parm01.addParm('nonHydrostatic',nonHydrostatic,PARM_BOOL);
  %%% exact volume conservation
  parm01.addParm('exactConserv',true,PARM_BOOL);
  %%% C-V scheme for Coriolis term
  parm01.addParm('useCDscheme',false,PARM_BOOL);  

  %%% file IO stuff
  parm01.addParm('readBinaryPrec',64,PARM_INT);
  parm01.addParm('useSingleCpuIO',true,PARM_BOOL);  
  parm01.addParm('debugLevel',-1,PARM_INT);
%   parm01.addParm('debugLevel',5,PARM_INT);

  %%%%% Vertical advection
  parm01.addParm('vectorInvariantMomentum',true,PARM_BOOL);
%   parm01.addParm('momImplVertAdv',true,PARM_BOOL);
%   parm01.addParm('tempImplVertAdv',true,PARM_BOOL);
%   parm01.addParm('saltImplVertAdv',true,PARM_BOOL);


  %%% Wet-point method at boundaries - may improve boundary stability
  parm01.addParm('useJamartWetPoints',true,PARM_BOOL);
  parm01.addParm('useJamartMomAdv',true,PARM_BOOL);

  %%% PARM02
  parm02.addParm('useSRCGSolver',true,PARM_BOOL);  
  parm02.addParm('cg2dMaxIters',1000,PARM_INT);  
  parm02.addParm('cg2dTargetResidual',1e-12,PARM_REAL);
  parm02.addParm('cg3dMaxIters',300,PARM_INT);  
  parm02.addParm('cg3dTargetResidual',1e-7,PARM_REAL);  
 
  %%% PARM03
  %parm03.addParm('alph_AB',1/2,PARM_REAL);
  %parm03.addParm('beta_AB',5/12,PARM_REAL);
  parm03.addParm('momDissip_In_AB',false,PARM_BOOL);
  
  %%%%%%% Added to make SEAICE work
  %parm03.addParm('momForcingOutAB',1,PARM_INT);
  parm03.addParm('tracForcingOutAB',1,PARM_INT);


  %%%%%%%
  parm03.addParm('nIter0',nIter0,PARM_INT);
  parm03.addParm('abEps',0.1,PARM_REAL);
  parm03.addParm('chkptFreq',1*t1day,PARM_REAL);
  parm03.addParm('pChkptFreq',5*t1day,PARM_REAL);
  parm03.addParm('taveFreq',0,PARM_REAL);
  parm03.addParm('dumpFreq',0*t1year,PARM_REAL);
  parm03.addParm('monitorFreq',1*t1day,PARM_REAL);
  parm03.addParm('cAdjFreq',0,PARM_REAL);
  parm03.addParm('dumpInitAndLast',true,PARM_BOOL);
  parm03.addParm('pickupStrictlyMatch',false,PARM_BOOL);
  if (periodicExternalForcing)
    parm03.addParm('periodicExternalForcing',periodicExternalForcing,PARM_BOOL);
    parm03.addParm('externForcingPeriod',externForcingPeriod,PARM_REAL);
    parm03.addParm('externForcingCycle',externForcingCycle,PARM_REAL);
  end
  
  %%% PARM04
  parm04.addParm('usingCartesianGrid',true,PARM_BOOL);
  parm04.addParm('usingSphericalPolarGrid',false,PARM_BOOL);    
 
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% GRID SPACING %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%    
    

  %%% Uniform meridional grid   
  dy = (Ly/Ny)*ones(1,Ny);
  yy = cumsum(dy);
  yy = yy-mean(yy);
  
  %%% Zonal grid
  if (use_3D)
    dx = Lx/Nx*ones(1,Nx);  
    xx = cumsum(dx);
    xx = xx-mean(xx);
  else
    dx = dy(1);
    xx = 0;
  end
 
  %%% Plotting mesh
  [Y,X] = meshgrid(yy,xx);
  
  %%% Grid spacing increases with depth, but spacings exactly sum to H
  zidx = 1:Nr;
  dz = H/Nr*ones(1,Nr);
  zz = -cumsum((dz+[0 dz(1:end-1)])/2);

  %%% Store grid spacings
  parm04.addParm('delX',dx,PARM_REALS);
  parm04.addParm('delY',dy,PARM_REALS);
  parm04.addParm('delR',dz,PARM_REALS);          
  
  
  
  
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%
  %%%%% BATHYMETRY %%%%%
  %%%%%%%%%%%%%%%%%%%%%%
   

  %%% Flat bottom  
  h = -H*ones(Nx,Ny);
  
  %%% Save as a parameter
  writeDataset(h,fullfile(inputpath,'bathyFile.bin'),ieee,prec);
  parm05.addParm('bathyFile','bathyFile.bin',PARM_STR); 
  
  
  
  
  
  
  
  

   
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% REFERENCE TEMPERATURE/SALINITY PROFILES %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
  %%% Quasi-tanh-shaped T/S profiles
  Zpyc = -150;
  Wpyc = 50;
  Smin = 34.3;
  Smax = 34.7;
  Tmin = 0.0901-0.0575*Smin; %%% MITgcm surface freezing temperature
  Tmax= 1;  
  gam_h = 0.01;
  tRef = Tmin + (Tmax-Tmin)*0.5*(1+0.5*(sqrt((1-(zz-Zpyc)/Wpyc).^2 + 4*gam_h*((zz-Zpyc)/Wpyc).^2)-sqrt((1+(zz-Zpyc)/Wpyc).^2 + 4*gam_h*((zz-Zpyc)/Wpyc).^2))/(1+gam_h)^(1/2));
  sRef = Smin + (Smax-Smin)*0.5*(1+0.5*(sqrt((1-(zz-Zpyc)/Wpyc).^2 + 4*gam_h*((zz-Zpyc)/Wpyc).^2)-sqrt((1+(zz-Zpyc)/Wpyc).^2 + 4*gam_h*((zz-Zpyc)/Wpyc).^2))/(1+gam_h)^(1/2)); 
  
  %%% Quasi-linear near-surface stratification
%   Zpyc = 0;
%   Hpyc = 300;  
%   Smin = 34.3;
%   Smax = 34.7;
%   Tmin = 0.0901-0.0575*Smin; %%% MITgcm surface freezing temperature
%   Tmax= 1;  
%   gam_h = 0.01;
%   tRef = Tmin + (Tmax-Tmin)*0.5*(sqrt((1-zz/Hpyc).^2 + 4*gam_h*(zz/Hpyc).^2)-sqrt((1+zz/Hpyc).^2 + 4*gam_h*(zz/Hpyc).^2))/(1+gam_h)^(1/2);
%   sRef = Smin + (Smax-Smin)*0.5*(sqrt((1-zz/Hpyc).^2 + 4*gam_h*(zz/Hpyc).^2)-sqrt((1+zz/Hpyc).^2 + 4*gam_h*(zz/Hpyc).^2))/(1+gam_h)^(1/2);
  
  %%% Plot the reference temperature
  if (showplots)
    figure(fignum);
    fignum = fignum + 1;
    clf;
    plot(tRef,zz);
    xlabel('\theta_r_e_f');
    ylabel('z','Rotation',0);
    title('Reference temperature');
  end
  
  %%% Plot the reference salinity
  if (showplots)
    figure(fignum);
    fignum = fignum + 1;
    clf;
    plot(sRef,zz);
    xlabel('S_r_e_f');
    ylabel('z','Rotation',0);
    title('Reference salinity');
  end
  
  
  
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DEFORMATION RADIUS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
  %%% Check Brunt-Vaisala frequency using full EOS
  pp = - zz;
%   SA = gsw_SA_from_SP(sRef,pp,-50,-64);  
%   CT = gsw_CT_from_pt(sRef,tRef);
%   [N2 pp_mid] = gsw_Nsquared(SA,CT,pp);
  [N2 pp_mid] = gsw_Nsquared(sRef,tRef,pp);
  dzData = zz(1:end-1)-zz(2:end);

  %%% Calculate internal wave speed and first Rossby radius of deformation
  N = sqrt(max(N2,0)); 
  Cig = sum(N.*dzData);
  Rd = Cig./(pi*abs(f0));

  if (showplots)
    figure(fignum);
    fignum = fignum + 1;
    clf;
    semilogx(N2,-pp_mid/1000);
    xlabel('N^2 (km)');
    ylabel('z (km)','Rotation',0);
    title('Buoyancy frequency');
  end
  
  if (showplots)
    ['First baroclinic Rossby deformation radius: ' num2str(Rd)]
  end
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% CALCULATE TIME STEP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
  
  
  %%% These estimates are in no way complete, but they give at least some
  %%% idea of the time step needed to keep things stable. In complicated 
  %%% simulations, preliminary tests may be required to estimate the
  %%% parameters used to calculate these time steps.        
  
  %%% Gravity wave CFL

  %%% Upper bound for absolute horizontal fluid velocity (m/s)
  %%% At the moment this is just an estimate
  Umax = 0.5;  
  Wmax = 0.1; %%% NEEDS TO BE INCREASED FOR DEEP PYCNOCLINES
  %%% Max gravity wave speed 
  cmax = max(Cig);
  %%% Max gravity wave speed using total ocean depth
  cgmax = Umax + cmax;
  %%% Gravity wave CFL
  deltaT_gw = min([0.5*dx/cmax,0.5*dy/cmax]);
  %%% Advective CFL
  deltaT_adv= min([0.5*dx/Umax,0.5*dy/Umax]);
  %%% Vertical advective CFL
  deltaT_vadv = min(0.5*dz/Wmax);
  %%% CFL time step based on full gravity wave speed
  deltaT_fgw = min([0.5*dx/cgmax,0.5*dy/cgmax]);
    
  %%% Other stability conditions
  
  %%% Inertial CFL time step (Sf0<=0.5)
  deltaT_itl = 0.5/abs(f0);
  %%% Time step constraint based on horizontal diffusion 
  deltaT_Ah = 0.5*min([dx dy])^2/(4*viscAh);    
  %%% Time step constraint based on vertical diffusion
  deltaT_Ar = 0.5*min(dz)^2 / (4*viscAr);  
  %%% Time step constraint based on biharmonic viscosity 
  deltaT_A4 = 0.5*min([dx dy])^4/(32*viscA4);
  %%% Time step constraint based on horizontal diffusion of temp 
  deltaT_KhT = 0.4*min([dx dy])^2/(4*diffKhT);    
  %%% Time step constraint based on vertical diffusion of temp 
  deltaT_KrT = 0.4*min(dz)^2 / (4*diffKrT);
  
  %%% Time step size  
  deltaT = min([deltaT_fgw deltaT_gw deltaT_adv deltaT_itl deltaT_Ah deltaT_Ar deltaT_KhT deltaT_KrT deltaT_A4]);
  if (nonHydrostatic)
    deltaT = min([deltaT deltaT_vadv]);
  end
  deltaT = round(deltaT);
  nTimeSteps = ceil(simTime/deltaT);
  simTimeAct = nTimeSteps*deltaT;
  
  %%% Write end time time step size  
  parm03.addParm('endTime',nIter0*deltaT+simTimeAct,PARM_INT);
  parm03.addParm('deltaT',deltaT,PARM_REAL); 

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% INITIAL DATA %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Random noise amplitude
  tNoise = 0.01;  
  sNoise = 0.001;
      
  %%% Align initial temp with background
  hydroTh = ones(Nx,Ny,Nr);
  hydroSa = ones(Nx,Ny,Nr);
  for k=1:1:Nr
    hydroTh(:,:,k) = squeeze(hydroTh(:,:,k))*tRef(k);
    hydroSa(:,:,k) = squeeze(hydroSa(:,:,k))*sRef(k);
  end
  
  %%% Add some random noise
  hydroTh = hydroTh + tNoise*(2*rand(Nx,Ny,Nr)-1);
  hydroSa = hydroSa + sNoise*(2*rand(Nx,Ny,Nr)-1);
  
  %%% Write to data files
  writeDataset(hydroTh,fullfile(inputpath,'hydrogThetaFile.bin'),ieee,prec); 
  parm05.addParm('hydrogThetaFile','hydrogThetaFile.bin',PARM_STR);
  writeDataset(hydroSa,fullfile(inputpath,'hydrogSaltFile.bin'),ieee,prec); 
  parm05.addParm('hydrogSaltFile','hydrogSaltFile.bin',PARM_STR);
  
  
  
  
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% TRACER DIFFUSION %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  diffK4T = diffK4Tgrid * max([dx dy])^4 / (32*deltaT);
  diffK4S = diffK4T;
  parm01.addParm('diffK4T',diffK4T,PARM_REAL); 
  parm01.addParm('diffK4S',diffK4S,PARM_REAL); 
 
  
  
  
    
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% SURFACE HEAT/SALT FLUXES %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  
  %%% Fixed salt flux in the BW formation region  
  salt_flux = zeros(Nx,Ny,nForcingPeriods);    
  tt = zeros(1,nForcingPeriods);
  for n=1:nForcingPeriods   
    t = (n-1)*externForcingPeriod;
    tt(n) = t;
    for j=1:Ny      
      for i=1:Nx             
        if (abs(yy(j)-Ypoly)<Wpoly/2)
          
          %%% Linear decrease of salt flux from salt_amp to zero after Tpoly days
%           salt_amp = -1e-1;
%           Tpoly = 3*t1day;
%           salt_flux(i,j,n) = salt_amp - salt_amp*min(1,t/Tpoly);             
          
          %%% Empirical ice growth equation from Anderson (1969)
          C1 = 1e4; %%% Constants from Anderson (1969) 
          C2 = 510;
          C3 = 6.7/86400;
          C4 = C2/(2*C1);
          C5 = C3/C1;
          Theta = 40; %%% Difference between atmospheric temperature and ocean freezing temperature 
          rho_i = 920; %%% Density of ice
          S_i = 5; %%% Salinity of ice
          e = -C4 + sqrt(C4^2+C5*t*Theta); %%% Ice thickness (m) as a function of time
          dedt = C3*Theta/(2*C1*e+C2); %%% Ice growth rate (m/s) as a function of time
          salt_flux(i,j,n) = -dedt*rho_i*(Smin-S_i); %%% Equivalent salt flux (g/m^2/s) as a function of time
          
        end
      end
    end         
  end  
  
 
  
  %%% Plot the surface salt flux
  if (showplots)
    figure(fignum);
    fignum = fignum + 1;
    clf;
%     contourf(X,Y,squeeze(salt_flux(:,:,1)));
    colorbar;
%     plot(yy,squeeze(salt_flux(1,:,1)));    
%     xlabel('y');
    plot(tt,squeeze(salt_flux(1,Ny/2,:)));    
    xlabel('t');
    ylabel('salt flux');
    title('Surface salt flux');
  end  
  
  %%% Save as parameters  
  writeDataset(salt_flux,fullfile(inputpath,'saltFluxFile.bin'),ieee,prec);
  parm05.addParm('saltFluxFile','saltFluxFile.bin',PARM_STR);
  












  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WRITE THE 'data' FILE %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Creates the 'data' file
  write_data(inputpath,PARM,listterm,realfmt);
  
 
  
  
  
  
  
  
  
  
  
  
  

  
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%
  %%%%% RBCS %%%%%
  %%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%
  
  
  
    
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%% RBCS SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% To store parameter names and values
  rbcs_parm01 = parmlist;
  rbcs_parm02 = parmlist;
  RBCS_PARM = {rbcs_parm01,rbcs_parm02};
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% RELAXATION PARAMETERS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  useRBCtemp = ~use_seaIce;
  useRBCsalt = false;
  useRBCuVel = false;
  useRBCvVel = false;
  tauRelaxT = 0.1*t1day;
  tauRelaxS = t1day;
  tauRelaxU = t1day;
  tauRelaxV = t1day;
  rbcs_parm01.addParm('useRBCtemp',useRBCtemp,PARM_BOOL);
  rbcs_parm01.addParm('useRBCsalt',useRBCsalt,PARM_BOOL);
  rbcs_parm01.addParm('useRBCuVel',useRBCuVel,PARM_BOOL);
  rbcs_parm01.addParm('useRBCvVel',useRBCvVel,PARM_BOOL);
  rbcs_parm01.addParm('tauRelaxT',tauRelaxT,PARM_REAL);
  rbcs_parm01.addParm('tauRelaxS',tauRelaxS,PARM_REAL);
  rbcs_parm01.addParm('tauRelaxU',tauRelaxU,PARM_REAL);
  rbcs_parm01.addParm('tauRelaxV',tauRelaxV,PARM_REAL);
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% RELAXATION TEMPERATURE %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  %%% Set relaxation temp equal to surface freezing temperature
  temp_relax = Tmin*ones(Nx,Ny,Nr);  

  %%% Save as parameters
  writeDataset(temp_relax,fullfile(inputpath,'sponge_temp.bin'),ieee,prec); 
  rbcs_parm01.addParm('relaxTFile','sponge_temp.bin',PARM_STR);  
  
  
  %%%%%%%%%%%%%%%%%%%%%  
  %%%%% RBCS MASK %%%%%
  %%%%%%%%%%%%%%%%%%%%%  
  
  
  %%% Mask is zero everywhere by default, i.e. no relaxation
  msk=zeros(Nx,Ny,Nr);  
    
  %%% Mask only nonzero at surface in the polynya
  for j=1:Ny      
    for i=1:Nx             
      if (abs(yy(j)-Ypoly)<Wpoly/2)
        msk(:,j,1) = 1;
      end
    end
  end         
  
  %%% Save as an input parameter
  writeDataset(msk,fullfile(inputpath,'rbcs_temp_mask.bin'),ieee,prec); 
  rbcs_parm01.addParm('relaxMaskFile(1)','rbcs_temp_mask.bin',PARM_STR);  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WRITE THE 'data.rbcs' FILE %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
 
  %%% Creates the 'data.rbcs' file
  write_data_rbcs(inputpath,RBCS_PARM,listterm,realfmt);
  
  
  
  
  
% %   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   %%%%%%%%% SEA ICE $ %%%%%%%%%%%
% %   
% %   
% %   % to store parameter names and values
%   seaice_parm01 = parmlist;
%   SEAICE_PARM = {seaice_parm01};
% % %   
% % %   
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%% SEA ICE  %%%%%%%%%%%%
%     %%%%%%%% PARAMETERS %%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Oringinal albedos from llc_1/48th  other values from llc_2160
% % 
%   SEAICEwriteState   = true;
%   SEAICEuseDYNAMICS  = true;
%   SEAICE_multDim     = 7;
% %   SEAICE_dryIceAlb   = 0.8783;
%   SEAICE_dryIceAlb   = 0.8509;
% %   SEAICE_wetIceAlb   = 0.7869;
%   SEAICE_wetIceAlb   = 0.7284;
% %   SEAICE_drySnowAlb  = 0.9482;
%   SEAICE_drySnowAlb  = 0.7754;
% %   SEAICE_wetSnowAlb  = 0.8216;
%   SEAICE_wetSnowAlb  = 0.7753;
%   SEAICE_waterDrag   = 5.5399;
%   SEAICE_drag        = 0.002;
% %   HO                 = 0.1;
%   HO                 = .05;
% 
% 
%   SEAICE_no_slip          = false;
% %   SEAICE_no_slip          = true;
% 
%   SEAICEadvScheme         = 7;
% %   SEAICEadvScheme         = 33;
% 
% 
%   %%%SOSEdoesn't have a seaice dataset for salinity, they used this value
%   %%%in their estimate
%   
%   LSR_ERROR               = 1.0e-4;
% %   LSR_ERROR               = 2.0e-4; 
%   MIN_ATEMP               = -40;
%   MIN_TICE                = -40;
%   SEAICE_area_reg         = 0.15;
%   SEAICE_hice_reg         = 0.1;
%   IMAX_TICE               = 6;
%   SEAICE_EPS		      = 1.0e-8;
% %   SEAICE_EPS              = 2.0e-9;
%   SEAICE_doOpenWaterMelt  = true;
%   SEAICE_areaLossFormula  = 1;
%   SEAICE_wetAlbTemp       = 0.0;
%   SEAICE_saltFrac         = 0.3;
% %   SEAICE_frazilFrac       = 0.003;
%  SEAICE_frazilFrac       = 0.01;
% %   SEAICE_frazilFrac       = 1.0;
% 
%   
%   seaice_parm01.addParm('LSR_ERROR',LSR_ERROR,PARM_REAL);
%   seaice_parm01.addParm('SEAICEwriteState',SEAICEwriteState,PARM_BOOL);
%   seaice_parm01.addParm('SEAICEuseDYNAMICS',SEAICEuseDYNAMICS,PARM_BOOL);
%   seaice_parm01.addParm('SEAICE_multDim',SEAICE_multDim,PARM_INT);
%   seaice_parm01.addParm('SEAICE_dryIceAlb',SEAICE_dryIceAlb,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_wetIceAlb',SEAICE_wetIceAlb,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_drySnowAlb',SEAICE_drySnowAlb,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_wetSnowAlb',SEAICE_wetSnowAlb,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_waterDrag',SEAICE_waterDrag,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_drag',SEAICE_drag,PARM_REAL);
%   seaice_parm01.addParm('HO',HO,PARM_REAL);
% %   seaice_parm01.addParm('SEAICE_dryIceAlb_south',SEAICE_dryIceAlb_south,PARM_REAL);
% %   seaice_parm01.addParm('SEAICE_wetIceAlb_south',SEAICE_wetIceAlb_south,PARM_REAL);
% %   seaice_parm01.addParm('SEAICE_drySnowAlb_south',SEAICE_drySnowAlb_south,PARM_REAL);
% %   seaice_parm01.addParm('SEAICE_wetSnowAlb_south',SEAICE_wetSnowAlb_south,PARM_REAL);
% %   seaice_parm01.addParm('SEAICE_waterDrag_south',SEAICE_waterDrag_south,PARM_REAL);
% %   seaice_parm01.addParm('SEAICE_drag_south',SEAICE_drag_south,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_no_slip',SEAICE_no_slip,PARM_BOOL);
% %   seaice_parm01.addParm('SEAICE_salinity',SEAICE_salinity,PARM_REAL);
%   seaice_parm01.addParm('SEAICEadvScheme',SEAICEadvScheme,PARM_INT);
%   seaice_parm01.addParm('MIN_ATEMP',MIN_ATEMP,PARM_REAL);
%   seaice_parm01.addParm('MIN_TICE',MIN_TICE,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_area_reg',SEAICE_area_reg,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_hice_reg',SEAICE_hice_reg,PARM_REAL);
%   seaice_parm01.addParm('IMAX_TICE',IMAX_TICE,PARM_INT);
%   seaice_parm01.addParm('SEAICE_EPS',SEAICE_EPS,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_doOpenWaterMelt',SEAICE_doOpenWaterMelt,PARM_BOOL);
%   seaice_parm01.addParm('SEAICE_areaLossFormula',SEAICE_areaLossFormula,PARM_INT);
%   seaice_parm01.addParm('SEAICE_wetAlbTemp',SEAICE_wetAlbTemp,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_saltFrac',SEAICE_saltFrac,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_frazilFrac',SEAICE_frazilFrac,PARM_REAL);
% 
%   
%   seaice_parm01.addParm('HeffFile',HeffFile,PARM_STR);
%   seaice_parm01.addParm('AreaFile',AreaFile,PARM_STR);
%   seaice_parm01.addParm('HsnowFile',HsnowFile,PARM_STR);
%   seaice_parm01.addParm('HsaltFile',HsaltFile,PARM_STR);
%   seaice_parm01.addParm('uIceFile',uIceFile,PARM_STR);
%   seaice_parm01.addParm('vIceFile',vIceFile,PARM_STR);
%   
%   
%      
% %   
% %   
% %   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   %%%%% WRITE THE 'data.seaice' FILE %%%
% %   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   
% %   
%   write_data_seaice(inputpath,SEAICE_PARM,listterm,realfmt);  
%   
%   
%   
%   
%   
%   
%   
%   
%   
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
%     %%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%
%     %%%%%%EXF PKG%%%%%
%     %%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%
%     
%     
%   %%% To store parameter names and values
%     
%      EXF_NML_01 = parmlist;
%      EXF_NML_02 = parmlist;
%      EXF_NML_03 = parmlist;
%      EXF_NML_04 = parmlist;
%      EXF_NML_OBCS = parmlist;
% 
%      EXF_PARM = {EXF_NML_01,EXF_NML_02,EXF_NML_03,EXF_NML_04,EXF_NML_OBCS};  
%     
%     
% 
% 
% % %     EXF_NML_01
% 
%   	exf_albedo        = 0.15;
%  	exf_scal_BulkCdn  = 1.015;
%  	exf_iprec         = 64;  
%  	useExfYearlyFields= false;
%  	useExfCheckRange  = false;
%  	useRelativeWind   = true;
% %  	useRelativeWind   = false;
%     repeatPeriod      = 283996800.0;
% %     exf_offset_atemp =  273.16;
%     
%     
% %%%runoff from ERA is in hours, need to convert to seconds
% %     exf_inscal_runoff = 1.14e-04;
%     
%     
%     apressurefile     = 'apressurefile.bin';
%     atempfile         = 'atempfile.bin';
%     aqhfile           = 'aqhfile.bin';
%     uwindfile         = 'uwindfile.bin';
%     vwindfile         = 'vwindfile.bin';
%     precipfile        = 'precipfile.bin';
%     swdownfile        = 'swdownfile.bin';
%     lwdownfile        = 'lwdownfile.bin';
% %     runofffile        = 'runofffile.bin';
%     
%     
% % "*period=-12" specifies monthly-mean forcing
%    apressurestartdate1 = 20070101;
%    apressurestartdate2 = 000000;
% %    apressureperiod     = 86400.0;
%    apressureperiod     = 10800.0;
%     
%     aqhstartdate1 = 20070101;
%     aqhstartdate2 = 000000;
% %     aqhperiod = 86400.0;
%     aqhperiod           = 10800.0;
% 
%     atempstartdate1 = 20070101;
%     atempstartdate2 = 000000;
% %     atempperiod = 86400.0;
%     atempperiod         = 10800.0;
%  
% 
%     uwindstartdate1 = 20070101;
%     uwindstartdate2 = 000000;
% %    uwindperiod = 86400.0;
%     uwindperiod         = 10800.0;
%  
%     vwindstartdate1 = 20070101;
%     vwindstartdate2 = 000000;
% %     vwindperiod = 86400.0;
%     vwindperiod         = 10800.0; 
%  
%     precipstartdate1 = 20070101;
%     precipstartdate2 = 000000;
% %     precipperiod = 86400.0;
%     precipperiod        = 10800.0; 
%  
% 
%     swdownstartdate1 = 20070101;
%     swdownstartdate2 = 000000;
% %     swdownperiod = 86400.0;
%     swdownperiod        = 10800.0;
% % 
% 
%     lwdownstartdate1 = 20070101;
%     lwdownstartdate2 = 000000;
% %     lwdownperiod = 86400.0;
%     lwdownperiod        = 10800.0;
% 
% 
% %     runoffstartdate1 = 20070101;
% %     runoffstartdate2 = 000000;
% %     runoffperiod = 2592000.0;
% %    runoffperiod        = 10800.0;
%  
%    precip_lon0 = xmin;
%    precip_lon_inc = dmxg(1);
%    precip_lat0 = ymin;
%    precip_lat_inc = dmyg;
%    precip_nlon = Nx;
%    precip_nlat = Ny;
%    
%    atemp_lon0 = xmin;
%    atemp_lon_inc = dmxg(1);
%    atemp_lat0 = ymin;
%    atemp_lat_inc = dmyg;
%    atemp_nlon = Nx;
%    atemp_nlat = Ny;
%    
%    apressure_lon0 = xmin;
%    apressure_lon_inc = dmxg(1);
%    apressure_lat0 = ymin;
%    apressure_lat_inc = dmyg;
%    apressure_nlon = Nx;
%    apressure_nlat = Ny;
%     
%    aqh_lon0 = xmin;
%    aqh_lon_inc = dmxg(1);
%    aqh_lat0 = ymin;
%    aqh_lat_inc = dmyg;
%    aqh_nlon = Nx;
%    aqh_nlat = Ny;
%    
%    uwind_lon0 = xmin;
%    uwind_lon_inc = dmxg(1);
%    uwind_lat0 = ymin;
%    uwind_lat_inc = dmyg;
%    uwind_nlon = Nx;
%    uwind_nlat = Ny;
%    
%    vwind_lon0 = xmin;
%    vwind_lon_inc = dmxg(1);
%    vwind_lat0 = ymin;
%    vwind_lat_inc = dmyg;
%    vwind_nlon = Nx;
%    vwind_nlat = Ny;
%    
%    swdown_lon0 = xmin;
%    swdown_lon_inc = dmxg(1);
%    swdown_lat0 = ymin;
%    swdown_lat_inc = dmyg;
%    swdown_nlon = Nx;
%    swdown_nlat = Ny;
%    
%    lwdown_lon0 = xmin;
%    lwdown_lon_inc = dmxg(1);
%    lwdown_lat0 = ymin;
%    lwdown_lat_inc = dmyg;
%    lwdown_nlon = Nx;
%    lwdown_nlat = Ny;
% 
% %    runoff_lon0 = xmin;
% %    runoff_lon_inc = dmxg(1);
% %    runoff_lat0 = ymin;
% %    runoff_lat_inc = dmyg;
% %    runoff_nlon = Nx;
% %    runoff_nlat = Ny;
%   
%  

%      
% %%%%%%%%%%%%%%%%%%%% ADD PARAMETERS
% 
%   EXF_NML_01.addParm('exf_albedo',exf_albedo,PARM_INT);
%   EXF_NML_01.addParm('exf_scal_BulkCdn',exf_scal_BulkCdn,PARM_REAL);
%   EXF_NML_01.addParm('exf_iprec',exf_iprec,PARM_INT);
%   EXF_NML_01.addParm('useExfYearlyFields',useExfYearlyFields,PARM_BOOL);
%   EXF_NML_01.addParm('useExfCheckRange',useExfCheckRange,PARM_BOOL);
%   EXF_NML_01.addParm('useRelativeWind',useRelativeWind,PARM_BOOL);
%   EXF_NML_01.addParm('repeatPeriod',repeatPeriod,PARM_REAL);
% %   EXF_NML_03.addParm('exf_offset_atemp',exf_offset_atemp,PARM_REAL);
% %   EXF_NML_03.addParm('exf_inscal_runoff',exf_inscal_runoff,PARM_REAL);
%   
%   
%   EXF_NML_02.addParm('apressurefile',apressurefile,PARM_STR);
%   EXF_NML_02.addParm('atempfile',atempfile,PARM_STR);
%   EXF_NML_02.addParm('aqhfile',aqhfile,PARM_STR);
%   EXF_NML_02.addParm('uwindfile',uwindfile,PARM_STR);
%   EXF_NML_02.addParm('vwindfile',vwindfile,PARM_STR);
%   EXF_NML_02.addParm('precipfile',precipfile,PARM_STR);
%   EXF_NML_02.addParm('swdownfile',swdownfile,PARM_STR);
%   EXF_NML_02.addParm('lwdownfile',lwdownfile,PARM_STR);
% %   EXF_NML_02.addParm('runofffile',runofffile,PARM_STR);
% 
% 
%   EXF_NML_02.addParm('apressurestartdate1',apressurestartdate1,PARM_INT);
%   EXF_NML_02.addParm('apressurestartdate2',apressurestartdate2,PARM_INT);
%   EXF_NML_02.addParm('apressureperiod',apressureperiod,PARM_REAL);
% 
%   EXF_NML_02.addParm('atempstartdate1',atempstartdate1,PARM_INT);
%   EXF_NML_02.addParm('atempstartdate2',atempstartdate2,PARM_INT);
%   EXF_NML_02.addParm('atempperiod',atempperiod,PARM_REAL);
%   
%   EXF_NML_02.addParm('aqhstartdate1',aqhstartdate1,PARM_INT);
%   EXF_NML_02.addParm('aqhstartdate2',aqhstartdate2,PARM_INT);
%   EXF_NML_02.addParm('aqhperiod',aqhperiod,PARM_REAL);
%   
%   EXF_NML_02.addParm('uwindstartdate1',uwindstartdate1,PARM_INT);
%   EXF_NML_02.addParm('uwindstartdate2',uwindstartdate2,PARM_INT);
%   EXF_NML_02.addParm('uwindperiod',uwindperiod,PARM_REAL);
%   
%   EXF_NML_02.addParm('vwindperiod',vwindperiod,PARM_REAL);
%   EXF_NML_02.addParm('vwindstartdate1',vwindstartdate1,PARM_INT);
%   EXF_NML_02.addParm('vwindstartdate2',vwindstartdate2,PARM_INT);
%   
%   EXF_NML_02.addParm('precipperiod',precipperiod,PARM_REAL);
%   EXF_NML_02.addParm('precipstartdate1',precipstartdate1,PARM_INT);
%   EXF_NML_02.addParm('precipstartdate2',precipstartdate2,PARM_INT);
%   
%   EXF_NML_02.addParm('swdownperiod',swdownperiod,PARM_REAL);
%   EXF_NML_02.addParm('swdownstartdate1',swdownstartdate1,PARM_INT);
%   EXF_NML_02.addParm('swdownstartdate2',swdownstartdate2,PARM_INT);
%   
%   EXF_NML_02.addParm('lwdownperiod',lwdownperiod,PARM_REAL);
%   EXF_NML_02.addParm('lwdownstartdate1',lwdownstartdate1,PARM_INT);
%   EXF_NML_02.addParm('lwdownstartdate2',lwdownstartdate2,PARM_INT);
%   
% %   EXF_NML_02.addParm('runoffperiod',runoffperiod,PARM_REAL);
% 
%   EXF_NML_04.addParm('precip_lon0',precip_lon0,PARM_REAL);
%   EXF_NML_04.addParm('precip_lon_inc',precip_lon_inc,PARM_REAL);
%   EXF_NML_04.addParm('precip_lat0',precip_lat0,PARM_REAL);
%   EXF_NML_04.addParm('precip_lat_inc',precip_lat_inc,PARM_REALS);
%   EXF_NML_04.addParm('precip_nlon',precip_nlon,PARM_INT);
%   EXF_NML_04.addParm('precip_nlat',precip_nlat,PARM_INT);
% 
%   EXF_NML_04.addParm('atemp_lon0',atemp_lon0,PARM_REAL);
%   EXF_NML_04.addParm('atemp_lon_inc',atemp_lon_inc,PARM_REAL);
%   EXF_NML_04.addParm('atemp_lat0',atemp_lat0,PARM_REAL);
%   EXF_NML_04.addParm('atemp_lat_inc',atemp_lat_inc,PARM_REALS);
%   EXF_NML_04.addParm('atemp_nlon',atemp_nlon,PARM_INT);
%   EXF_NML_04.addParm('atemp_nlat',atemp_nlat,PARM_INT);
% 
%   EXF_NML_04.addParm('apressure_lon0',apressure_lon0,PARM_REAL);
%   EXF_NML_04.addParm('apressure_lon_inc',apressure_lon_inc,PARM_REAL);
%   EXF_NML_04.addParm('apressure_lat0',apressure_lat0,PARM_REAL);
%   EXF_NML_04.addParm('apressure_lat_inc',apressure_lat_inc,PARM_REALS);
%   EXF_NML_04.addParm('apressure_nlon',apressure_nlon,PARM_INT);
%   EXF_NML_04.addParm('apressure_nlat',apressure_nlat,PARM_INT);
% 
%   EXF_NML_04.addParm('aqh_lon0',aqh_lon0,PARM_REAL);
%   EXF_NML_04.addParm('aqh_lon_inc',aqh_lon_inc,PARM_REAL);
%   EXF_NML_04.addParm('aqh_lat0',aqh_lat0,PARM_REAL);
%   EXF_NML_04.addParm('aqh_lat_inc',aqh_lat_inc,PARM_REALS);
%   EXF_NML_04.addParm('aqh_nlon',aqh_nlon,PARM_INT);
%   EXF_NML_04.addParm('aqh_nlat',aqh_nlat,PARM_INT);
% 
%   EXF_NML_04.addParm('uwind_lon0',uwind_lon0,PARM_REAL);
%   EXF_NML_04.addParm('uwind_lon_inc',uwind_lon_inc,PARM_REAL);
%   EXF_NML_04.addParm('uwind_lat0',uwind_lat0,PARM_REAL);
%   EXF_NML_04.addParm('uwind_lat_inc',uwind_lat_inc,PARM_REALS);
%   EXF_NML_04.addParm('uwind_nlon',uwind_nlon,PARM_INT);
%   EXF_NML_04.addParm('uwind_nlat',uwind_nlat,PARM_INT);
% 
%   EXF_NML_04.addParm('vwind_lon0',vwind_lon0,PARM_REAL);
%   EXF_NML_04.addParm('vwind_lon_inc',vwind_lon_inc,PARM_REAL);
%   EXF_NML_04.addParm('vwind_lat0',vwind_lat0,PARM_REAL);
%   EXF_NML_04.addParm('vwind_lat_inc',vwind_lat_inc,PARM_REALS);
%   EXF_NML_04.addParm('vwind_nlon',vwind_nlon,PARM_INT);
%   EXF_NML_04.addParm('vwind_nlat',vwind_nlat,PARM_INT);
% 
%   EXF_NML_04.addParm('swdown_lon0',swdown_lon0,PARM_REAL);
%   EXF_NML_04.addParm('swdown_lon_inc',swdown_lon_inc,PARM_REAL);
%   EXF_NML_04.addParm('swdown_lat0',swdown_lat0,PARM_REAL);
%   EXF_NML_04.addParm('swdown_lat_inc',swdown_lat_inc,PARM_REALS);
%   EXF_NML_04.addParm('swdown_nlon',swdown_nlon,PARM_INT);
%   EXF_NML_04.addParm('swdown_nlat',swdown_nlat,PARM_INT);
% 
%   EXF_NML_04.addParm('lwdown_lon0',lwdown_lon0,PARM_REAL);
%   EXF_NML_04.addParm('lwdown_lon_inc',lwdown_lon_inc,PARM_REAL);
%   EXF_NML_04.addParm('lwdown_lat0',lwdown_lat0,PARM_REAL);
%   EXF_NML_04.addParm('lwdown_lat_inc',lwdown_lat_inc,PARM_REALS);
%   EXF_NML_04.addParm('lwdown_nlon',lwdown_nlon,PARM_INT);
%   EXF_NML_04.addParm('lwdown_nlat',lwdown_nlat,PARM_INT);
% 
% %   EXF_NML_04.addParm('runoff_lon0',runoff_lon0,PARM_REAL);
% %   EXF_NML_04.addParm('runoff_lon_inc',runoff_lon_inc,PARM_REAL);
% %   EXF_NML_04.addParm('runoff_lat0',runoff_lat0,PARM_REAL);
% %   EXF_NML_04.addParm('runoff_lat_inc',runoff_lat_inc,PARM_REALS);
% %   EXF_NML_04.addParm('runoff_nlon',runoff_nlon,PARM_INT);
% %   EXF_NML_04.addParm('runoff_nlat',runoff_nlat,PARM_INT);
% 
% 
%   %%z% Create the data.exf file
%   write_data_exf(inputpath,EXF_PARM,listterm,realfmt);






  
  
  
  
  
  
 
 
  
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DIAGNOSTICS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
    
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DIAGNOSTICS SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% To store parameter names and values
  diag_parm01 = parmlist;
  diag_parm02 = parmlist;
  DIAG_PARM = {diag_parm01,diag_parm02};
  diag_matlab_parm01 = parmlist;
  DIAG_MATLAB_PARM = {diag_matlab_parm01}; %%% Matlab parameters need to be different
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DIAGNOSTICS PARAMETERS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  
  %%% Stores total number of diagnostic quantities
  ndiags = 0;
  
  diag_fields_avg = ...
  {...

%     'UV_VEL_Z','WU_VEL','WV_VEL'... %%% Momentum fluxes
%     'UVELSLT','VVELSLT','WVELSLT', ... %%% Salt fluxes

%     'PHIHYD', ... %%% Pressure
%     'LaUH1TH','LaVH1TH','LaHw1TH','LaHs1TH' ... %%% LAYERS fluxes
%     'ADVr_TH','ADVx_TH','ADVy_TH', ... %%%%% Advective fluxes of H.B.
%     'KPPg_TH', ... %%%% KPP non-local flux of P.T.
%     'DFrE_TH','DFrI_TH','DFxE_TH','DFyE_TH', ... %%%% Diffusive Fluxes
%     'oceQsw','TOTTTEND','WTHMASS','TFLUX'...
%      'UVELSQ','VVELSQ','WVELSQ',... %%% For Kinetic Energy
    %%'ETAN',...
%    'SHIgammT','SHIgammS','SHIuStar','SHI_mass', ... %%%%% Ice shelf melt
%      'SIarea','SIheff','SIhsnow','SIhsalt','SIuice','SIvice'...%%% Sea ice state
%      'UVEL','VVEL','WVEL'... %%% Velocities
%      'THETA', ... %%% Temperature
%      'SALT', ... %%% Salinity
%      'SHIfwFlx'...%%'SHIhtFlx','SHIForcT','SHIForcS', ...
%      'UVELTH','VVELTH','WVELTH', ... %%% Temperature fluxes
  };
  
  numdiags_avg = length(diag_fields_avg);  
  diag_freq_avg = 5*t1day;
  diag_phase_avg = 0;    
     
  diag_parm01.addParm('diag_mnc',true,PARM_BOOL);  
  for n=1:numdiags_avg    
    
    ndiags = ndiags + 1;
    
    diag_parm01.addParm(['fields(1,',num2str(n),')'],diag_fields_avg{n},PARM_STR);  
    diag_parm01.addParm(['fileName(',num2str(n),')'],diag_fields_avg{n},PARM_STR);  
    diag_parm01.addParm(['frequency(',num2str(n),')'],diag_freq_avg,PARM_REAL);  
    diag_parm01.addParm(['timePhase(',num2str(n),')'],diag_phase_avg,PARM_REAL); 
    diag_matlab_parm01.addParm(['diag_fields{1,',num2str(n),'}'],diag_fields_avg{n},PARM_STR);  
    diag_matlab_parm01.addParm(['diag_fileNames{',num2str(n),'}'],diag_fields_avg{n},PARM_STR);  
    diag_matlab_parm01.addParm(['diag_frequency(',num2str(n),')'],diag_freq_avg,PARM_REAL);  
    diag_matlab_parm01.addParm(['diag_timePhase(',num2str(n),')'],diag_phase_avg,PARM_REAL);  
    
  end
  
  diag_fields_inst = ...
  {...  
    'ETAN', ...  %%% SSH
    'UVEL','VVEL','WVEL'...%%%velocities
    'THETA', ... %%% Temperature
    'SALT', ... %%% Salinity
     %'KPPhbl' ... %%% BL depth
    }; 
  
  numdiags_inst = length(diag_fields_inst);  
  if (use_3D)
    diag_freq_inst = 0.1*t1day;
  else
    diag_freq_inst = 0.01*t1day;
  end
  diag_phase_inst = 0;
  
  for n=1:numdiags_inst    
    
    ndiags = ndiags + 1;
    
    diag_parm01.addParm(['fields(1,',num2str(ndiags),')'],diag_fields_inst{n},PARM_STR);  
    diag_parm01.addParm(['fileName(',num2str(ndiags),')'],[diag_fields_inst{n},'_inst'],PARM_STR);  
    diag_parm01.addParm(['frequency(',num2str(ndiags),')'],-diag_freq_inst,PARM_REAL);  
    diag_parm01.addParm(['timePhase(',num2str(ndiags),')'],diag_phase_inst,PARM_REAL); 
    diag_matlab_parm01.addParm(['diag_fields(1,',num2str(ndiags),')'],diag_fields_inst{n},PARM_STR);  
    diag_matlab_parm01.addParm(['diag_fileNames(',num2str(ndiags),')'],[diag_fields_inst{n},'_inst'],PARM_STR);  
    diag_matlab_parm01.addParm(['diag_frequency(',num2str(ndiags),')'],-diag_freq_inst,PARM_REAL);  
    diag_matlab_parm01.addParm(['diag_timePhase(',num2str(ndiags),')'],diag_phase_inst,PARM_REAL);  
    
  end
  
  %%% Create the data.diagnostics file
  write_data_diagnostics(inputpath,DIAG_PARM,listterm,realfmt);
  
  %%% Create the DIAGNOSTICS_SIZE.h file
  createDIAGSIZEh(codepath,ndiags,Nr);
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%
  %%%%% PACKAGES %%%%%
  %%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%
  
  packages = parmlist;
  PACKAGE_PARM = {packages};  
  
  packages.addParm('useDiagnostics',true,PARM_BOOL);    
  packages.addParm('useKPP',~nonHydrostatic,PARM_BOOL);
  packages.addParm('useRBCS',~use_seaIce,PARM_BOOL);      
  packages.addParm('useSEAICE',use_seaIce,PARM_BOOL);  
  packages.addParm('useEXF',use_seaIce,PARM_BOOL);  
  
  %%% Create the data.pkg file
  write_data_pkg(inputpath,PACKAGE_PARM,listterm,realfmt);
  
  
  
  
  
  
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WRITE PARAMETERS TO A MATLAB FILE %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Creates a matlab file defining all input parameters
  ALL_PARMS = [PARM PACKAGE_PARM DIAG_MATLAB_PARM];
  if (~nonHydrostatic)
    ALL_PARMS = [ALL_PARMS KPP_PARM];
  end
  if (use_seaIce)
    ALL_PARMS = [ALL_PARMS SEAICE_PARM EXF_PARM];
  else
    ALL_PARMS = [ALL_PARMS RBCS_PARM];
  end    
  write_matlab_params(inputpath,ALL_PARMS,realfmt);
  
end
