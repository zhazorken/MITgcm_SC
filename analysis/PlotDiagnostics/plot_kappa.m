%%%
%%% plot time-averaged eddy diffusivity diagnostics
%%%

addpath ./data_files

R0 = 2*1000; % initial eddy radius in meters (approxmiated based on results from Cohanim et al. (2020))
vmax = 0.1*86400; % rough estimate of maximum horizontal velocity in m/day
t_break = 10; % time in days where eddy growth is approxmiately dominated by merging 

tau = 5*ones(1,48);

t = [1:48];
tau(1:48) = round(4*(R0+3/5*t*1000)/vmax);
kappa_avg = zeros(size(kappa));

for f=1:size(kappa,1)
    for n=1:length(tau)
        knum = movmean(squeeze(k_num(f,:,:)),tau(n),1);
        kdenom = movmean(squeeze(k_denom(f,:,:)),tau(n),1);
        kappa_avg(f,n,:) = knum(n,:)./kdenom(n,:);
    end
end

scrsz = get(0,'ScreenSize');
framepos = [0.3*1 0.3*1 1200 900]; %[0.25*scrsz(3) 0.15*scrsz(4) 900 400];
figure;
clf;
set(gcf,'Position',framepos);
set(gcf,'Color','w');

pcolor(abs(squeeze(kappa_avg(1,:,:)))')
shading interp;
cmocean('-deep',256);
colorbar; 
caxis([0 100]);
xlabel('Time (days)')
ylabel('y (km)')
ax = gca;
ax.FontSize = 28;
ax.FontName = 'Times New Roman';
ax.YGrid = 'on';
ax.Box = 'on';
ax.LineWidth = 2;