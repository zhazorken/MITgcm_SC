%%%
%%% function to produce a 3D snapshot

%function CreateSnapshot3D(xx,yy,zz,pdens,uvel,vvel,Nx,Ny,Nr,ff,delY,delX)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Vorticity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath ./cm_and_cb_utilities/cm_and_cb_utilities
addpath ./data_files

scrsz = get(0,'ScreenSize');
framepos = [0.3*1 0.3*1 1200 900]; %[0.25*scrsz(3) 0.15*scrsz(4) 900 400];
figure;
clf;
set(gcf,'Position',framepos);
set(gcf,'Color','w');
vid=VideoWriter('Ref2');
vid.FrameRate=3;
open(vid);

for n=1:48
    vort = squeeze(vort3D(n,:,:,:));
    pdens = squeeze(pdens3D(n,:,:,:));
    yy = xx;
    %%% Create 3D vorticity field
    [Y,X] = meshgrid(yy,xx);
%     vort = zeros(Nx,Ny,Nr);
%     vort(:,2:Ny,:) = - (uvel(:,2:Ny,:)-uvel(:,1:Ny-1,:))/delY(1);
%     vort = vort + (vvel([2:Nx 1],:,:)-vvel(:,:,:))/delX(1);
    ff = f0*ones(Nx,Ny);
    vort=vort./abs(ff);
    
    % Vorticity surface
    surf_lev=3;
    pvort=surf(X/1000,Y/1000,zz(surf_lev)*ones(size(X)),vort(:,:,surf_lev));
    pvort.FaceColor='texturemap';
    colormap(cmocean('balance',256));
    caxis(gca,[-1 1]);
    %caxis('manual')
    shading interp
    pvort.EdgeColor='none';
    colorbar;
    alpha(pvort,0.55);
    %freezeColors
    hold on;
    % % Vorticity Isosurface
    % [Y X Z]=meshgrid(yy(3:end),xx(3:end),zz(3:end));
    % isovort = vort(3:end,3:end,3:end);
    % %isopyc(:,:,1:round(Hm_lev/2)) = nan;
    % vort_plot=0.3; % isovort(0.1,0.1,:);
    % fv=isosurface(X/1000,Y/1000,Z,isovort,vort_plot);
    % p=patch(fv);
    % p.FaceColor=[208/255,49/255,45/255];
    % p.EdgeColor='none';
    % alpha(p,0.1);
    % freezeColors
    % hold on;
    %
    % [Y X Z]=meshgrid(yy(3:end),xx(3:end),zz(3:end));
    % isovort = vort(3:end,3:end,3:end);
    % vort_plot=-0.3; % isovort(0.1,0.1,:);
    % fv=isosurface(X/1000,Y/1000,Z,isovort,vort_plot);
    % p=patch(fv);
    % p.FaceColor=[58/255,67/255,186/255];
    % p.EdgeColor='none';
    % alpha(p,0.2);
    % freezeColors
    % hold on;
    
    % [Y X Z]=meshgrid(yy(2:end),xx(2:end),zz(Hm_lev-2:Hm_lev));
    % isovort = vort(2:end,2:end,Hm_lev-2:Hm_lev);
    % vort_plot=0.1;
    % fv=isosurface(X/1000,Y/1000,Z,isovort,vort_plot);
    % p=patch(fv);
    % p.FaceColor=[58/255,67/255,186/255];
    % p.EdgeColor='none';
    % alpha(p,0.2);
    % freezeColors
    % hold on;
    
    % [Y X Z]=meshgrid(yy(2:end),xx(2:end),zz(Hm_lev-2:Hm_lev));
    % isovort = vort(2:end,2:end,Hm_lev-2:Hm_lev);
    % vort_plot=-0.1;
    % fv=isosurface(X/1000,Y/1000,Z,isovort,vort_plot);
    % p=patch(fv);
    % p.FaceColor=[58/255,67/255,186/255];
    % p.EdgeColor='none';
    % alpha(p,0.2);
    % freezeColors
    % hold on;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Density %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set isopycnal depth
    Hm_lev=18;
    
    % Plot density at Hm isopycnal
    [Y X Z]=meshgrid(yy,xx,zz);
    isopyc = pdens;
    isopyc(:,:,1:surf_lev) = nan;
    pdens_plot=isopyc(200,200,Hm_lev);
    fv=isosurface(X/1000,Y/1000,Z,isopyc,pdens_plot);
    p=patch(fv);
    p.FaceColor=[79 66 181]/255;
    p.EdgeColor='none';
    alpha(p,0.75);
    freezeColors
    hold on;
    
    %%% Plot background
    pdens_end = squeeze(pdens(2,:,:));
    pdens_end(pdens_end==0) = NaN;
    [ZZZ,YYY] = meshgrid(zz,yy);
    zz_i = zz(end)+0.1/2:0.1:zz(1)-0.1/2;
    Nz_i = length(zz_i);
    [ZZZ_i,YYY_i] = meshgrid(zz_i,yy);
    pdens_i = zeros(Ny,Nz_i);
    for j=1:Ny
        pdens_i(j,:) = interp1(zz,pdens_end(j,:),zz_i,'linear');
    end
    p2 = surface(-24*ones(size(YYY_i)),YYY_i/1000,ZZZ_i,pdens_i,'FaceAlpha',0.9,'EdgeColor','none');
    p2.FaceColor = 'texturemap';
    colormap(cmocean('-haline',150));
    caxis([1027.4 1027.5])
    hold on;
    %contour(pdens_i,70);
    freezeColors;
    hold on;
    
    %%% Plot background
    pdens_end = squeeze(pdens(:,2,:));
    pdens_end(pdens_end==0) = NaN;
    [ZZZ,XXX] = meshgrid(zz,xx);
    zz_i = zz(end)+0.1/2:0.1:zz(1)-0.1/2;
    Nz_i = length(zz_i);
    [ZZZ_i,XXX_i] = meshgrid(zz_i,xx);
    pdens_i = zeros(Ny,Nz_i);
    for j=1:Nx
        pdens_i(j,:) = interp1(zz,pdens_end(j,:),zz_i,'linear');
    end
    %ax1 = axes;
    p2 = surface(XXX_i/1000,-24*ones(size(XXX_i)),ZZZ_i,pdens_i,'FaceAlpha',0.9,'EdgeColor','none');
    p2.FaceColor = 'texturemap';
    colormap(cmocean('-haline',150));
    %colormap(pmkmp(256,'Swtth'))
    caxis([1027.4 1027.5])
    freezeColors
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Decorations %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    az = 125;
    elev = 15;
    view(az,elev); %view(125,21) %view(-113,-19); %view(62,-25);
    axis tight;
    xlim([-24 25]);
    ylim([-24 25]);
    zlim([-120 -2.5])
    handle = gca;
    set(handle, 'fontsize',30)
    title(['t = ',num2str(n),' days'],'interpreter','latex','fontsize',48);
    zlabel('z (m)','interpreter','latex','fontsize',48);
    ylabel('y (km)','interpreter','latex','fontsize',48);
    xlabel('x (km)','interpreter','latex','fontsize',48);
    set(handle,'Position',[ 0.14   0.1100    0.7    0.8]);
    handle.Box = 'on';
    handle.LineWidth = 2;
    cb = colorbar;
    set(cb,'Position',[0.89    0.3    0.02    0.47],'YTick',[1027.4 1027.425 1027.45 1027.475 1027.5]); %[0.9199    0.6983    0.0201    0.2117]
    cb.TickLength = 0.025;
    cb.FontSize = 24;
    cb.Box = 'on';
    cb.LineWidth = 2.5;
    colorTitle = get(cb, 'Title');
    set(colorTitle,'String','$\sigma \mathrm{(kg \, m^{-3})}$','interpreter','latex','fontsize',28)
    %pbaspect([1.25 1.25 1]);
    %annotation('textbox',[0.73 0.9 0.15 0.01],'String',{'$\zeta/|f|$'},'fontsize',fontsize+2,'LineStyle','None','interpreter','latex');
    camlight('headlight');
    lightangle(az,45);
    %lightangle(47,22);
    lighting gouraud;
    %savefig('ref3D');
    
    hold off;
    
   M(n) = getframe(gcf);
    writeVideo(vid,M(n))
end

close(vid)
save('3Dvid','Ref.avi')