%%%
%%% Potential density field, vorticity and  buoyacy sensitivty plots
%%%

addpath ./mat_files/
scrsz = get(0,'ScreenSize');
framepos = [0.3*1 0.3*1 1200 900]; %[0.25*scrsz(3) 0.15*scrsz(4) 900 400];
figure;
clf;
set(gcf,'Position',framepos);
set(gcf,'Color','w');

subplot_rows = 4; % number of runs
subplot_cols = 3;

%%% Choose Runs
%f = [2 9 21 25]; % specify run number from Exps


sp_count = 0;
day = 10; % specify day to plot vorticity and potential density field

for i=1:subplot_rows    %%% Vorticity Field
    
    % plot vorticity field
    subplot(subplot_rows,subplot_cols,sp_count+3)
    [YY,XX] = meshgrid(yy,xx);
    if i==1
        vort = squeeze(Vort(i,:,:));
    else
        vort = squeeze(Vort(i+1,:,:));
    end
    pcolor(XX/1000,YY/1000,vort)
    shading interp;
    colormap(cmocean('balance',256));
    caxis([-1 1]*10^(-3));
    pbaspect([1 1 1])
    ax2 = gca;
    ax2.FontSize = 28;
    ax2.FontName = 'Times New Roman';
    ax2.Box = 'on';
    ylabel('y (km)');
    sp_count = sp_count+3;
end

sp_count = 0;
for i=1:subplot_rows
    %%% Potential Density Field
    [ZZ,XX] = meshgrid(zz,xx);
    subplot(subplot_rows,subplot_cols,sp_count+1)
    pdens = squeeze(Pdens(f(i),day,:,:));
    pcolor(XX/1000,ZZ,pdens); 
    shading interp; 
    cmocean('-haline',256);
    caxis([1027.4 1027.6])
    ylim([-350 -2.5]);
    pbaspect([1.75 1 1]);
    hold on;
    ax = gca;
    ax.FontSize = 28;
    ax.FontName = 'Times New Roman';
    ax.Box = 'on';
    ylabel('z (m)');
    sp_count = sp_count+3;
    
end

sp_count = 0;
 for i=1:subplot_rows   
    %%% Buoyancy
    subplot(subplot_rows,subplot_cols,sp_count+2)
    b = squeeze(B(f(i),:,:));
    pcolor(b');
    shading interp;
    cmocean('rain',256)
    caxis([-0.06 0])
    pbaspect([1.75 1 1]);
    ax3 = gca;
    ax3.FontSize = 28;
    ax3.FontName = 'Times New Roman';
    ax3.Box = 'on';
    sp_count = sp_count+3;
end