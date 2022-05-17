%%%
%%% avg_t.m
%%%
%%% Calculates the time average of the output fields from MITgcm runs.
%%%

%%% Read experiment data
clear diag_fields;
clear diag_timePhase;
clear diag_fileNames;
clear diag_frequency;
loadexp;

%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics.
dumpFreq = diag_frequency(1);
nDumps = floor(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);

%%% Calculate time average for whichever fields are present
for m=1:length(diag_fields)
  if (m == 1) %%% N.B. This won't work if the first field isn't one of those listed below
    flag = '';
  else    
    flag = '-append';
  end
  if (strcmp(diag_fields(m),'UVEL'));
    uu = readIters(exppath,'UVEL',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);    
    save([expname,'_tavg.mat'],'uu',flag);
  end
  if (strcmp(diag_fields(m),'VVEL'));
    vv = readIters(exppath,'VVEL',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
    save([expname,'_tavg.mat'],'vv',flag);
  end
  if (strcmp(diag_fields(m),'WVEL'));
    ww = readIters(exppath,'WVEL',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
    save([expname,'_tavg.mat'],'ww',flag);
  end
  if (strcmp(diag_fields(m),'THETA'));    
    tt = readIters(exppath,'THETA',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
    save([expname,'_tavg.mat'],'tt',flag);
  end
  if (strcmp(diag_fields(m),'SALT'));
    ss = readIters(exppath,'SALT',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
    save([expname,'_tavg.mat'],'ss',flag);
  end
  if (strcmp(diag_fields(m),'PHIHYD'));
    pp = readIters(exppath,'PHIHYD',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
    save([expname,'_tavg.mat'],'pp',flag);
  end
  if (strcmp(diag_fields(m),'UVELTH'));
    ut = readIters(exppath,'UVELTH',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
    save([expname,'_tavg.mat'],'ut',flag);
  end
  if (strcmp(diag_fields(m),'VVELTH'));
    vt = readIters(exppath,'VVELTH',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
    save([expname,'_tavg.mat'],'vt',flag);
  end
  if (strcmp(diag_fields(m),'WVELTH'));
    wt = readIters(exppath,'WVELTH',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
    save([expname,'_tavg.mat'],'wt',flag);
  end
  if (strcmp(diag_fields(m),'UVELSLT'));
    us = readIters(exppath,'UVELSLT',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
    save([expname,'_tavg.mat'],'us',flag);
  end
  if (strcmp(diag_fields(m),'VVELSLT'));
    vs = readIters(exppath,'VVELSLT',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
    save([expname,'_tavg.mat'],'vs',flag);
  end
  if (strcmp(diag_fields(m),'WVELSLT'));
    ws = readIters(exppath,'WVELSLT',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
    save([expname,'_tavg.mat'],'ws',flag);
  end
  if (strcmp(diag_fields(m),'THETASQ'));
    tsq = readIters(exppath,'THETASQ',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
    save([expname,'_tavg.mat'],'tsq',flag);
  end
  if (strcmp(diag_fields(m),'SALTSQ'));
    ssq = readIters(exppath,'SALTSQ',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
    save([expname,'_tavg.mat'],'ssq',flag);
  end
  if (strcmp(diag_fields(m),'THSLT'));
    ts = readIters(exppath,'THSLT',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
    save([expname,'_tavg.mat'],'ts',flag);
  end
  if (strcmp(diag_fields(m),'UVELSQ'));
    usq = readIters(exppath,'UVELSQ',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
    save([expname,'_tavg.mat'],'usq',flag);
  end
  if (strcmp(diag_fields(m),'VVELSQ'));
    vsq = readIters(exppath,'VVELSQ',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
    save([expname,'_tavg.mat'],'vsq',flag);
  end
  if (strcmp(diag_fields(m),'WVELSQ'));
    wsq = readIters(exppath,'WVELSQ',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
    save([expname,'_tavg.mat'],'wsq',flag);
  end
  if (strcmp(diag_fields(m),'UV_VEL_Z'));
    uv = readIters(exppath,'UV_VEL_Z',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
    save([expname,'_tavg.mat'],'uv',flag);
  end
  if (strcmp(diag_fields(m),'WU_VEL'));
    uw = readIters(exppath,'WU_VEL',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
    save([expname,'_tavg.mat'],'uw',flag);
  end
  if (strcmp(diag_fields(m),'WV_VEL'));
    vw = readIters(exppath,'WV_VEL',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
    save([expname,'_tavg.mat'],'vw',flag);
  end
end
