%%%
%%% Plot diagnostics of eddy size for reference case
%%%
%%% Subplot 1: plotting KE as a function of wavenumber for specific days. Scatters data and plots smoothed lines.
%%%
%%% Subplot 2: plot of R(t). Scatters data and plots smooth lines with
%%% markers corresponding to the days in subplot 1. 
%%%

scrsz = get(0,'ScreenSize');
framepos = [0.3*1 0.3*1 1200 900]; %[0.25*scrsz(3) 0.15*scrsz(4) 900 400];
figure;
clf;
set(gcf,'Position',framepos);
set(gcf,'Color','w');

%%% Subplot 1
subplot(1,2,1);
T = [3 5 10 48];
colors = zeros(length(T),3);
colors(1,:) = [0.00,0.45,0.74];
colors(2,:) = [0.49,0.18,0.56];
colors(3,:) = [0.47,0.67,0.19];
colors(4,:) = [0.93,0.69,0.13];
ccount=0;
for n=T
    ccount = ccount+1;
    color = colors(ccount,:);
    ke = log10(squeeze(KE(1,n,:)/sum(KE(1,n,:))));
    p = scatter([1:30],ke(1:30),'x');
    p.MarkerFaceColor = color;
    %p.MarkerEdgeColor = 1; 
    p.MarkerFaceAlpha= 1;
    p.SizeData = 70;
    p.LineWidth = 2;
    hold on;
    p = plot(smooth(ke(1:30),3),'linewidth',4,'color',color); % plot log10 of normalized KE
    p.Color(4) = 0.5;
    cent = Centroid(1,n);
    plot(ones(1,size(KE,3))*cent,linspace(-7,0,size(KE,3)),':','color',color,'linewidth',2.5)
end

pbaspect([1 1 1])
ylim([-7 0])
xlim([1 30])
xlabel('wavenumber,k','FontName','Times New Roman','interpreter','latex')
ylabel('$\mathrm{log_{10}}(\hat{E}$)','interpreter','latex')
ax = gca;
ax.FontSize = 28;
ax.YGrid = 'on';
ax.Box = 'on';
ax.LineWidth = 2;
%ax.YTickLabel = [10^-7 10^-6 10^-5 10^-4 10^-3 10^-2 10^-1 10^0];

%%% Subplot 2
subplot(1,2,2);
p = scatter([1:48],R(1,:));
p.MarkerFaceColor = colors(1,:);
p.MarkerEdgeColor = 'none';
p.MarkerFaceAlpha= 0.5;
hold on;
plot([1:48],smooth(R(1,:)),'color',colors(1,:),'linewidth',4); hold on;
for n=1:size(T,2)
    scatter(T(n),R(1,T(n)),'s','MarkerEdgeColor','none','MarkerFaceColor',[210/255,21/255,2/255]); hold on
end
pbaspect([1 1 1])
ylim([0 25])
xlim([1 48])
xlabel('Time (days)','interpreter','latex')
ylabel('Eddy Size (km)','interpreter','latex')
ax = gca;
ax.FontSize = 28;
ax.YGrid = 'on';
ax.XGrid = 'on';
ax.Box = 'on';
ax.LineWidth = 2;
