%%%
%%%
%%% Plot eddy size as a function fo time for all runs. Need eddys size .mat file.
%%%
%%%

T = 30; % choose cut-off time
expNum = size(R,1); % length of runs

scrsz = get(0,'ScreenSize');
framepos = [0.3*1 0.3*1 1200 900]; %[0.25*scrsz(3) 0.15*scrsz(4) 900 400];
figure;
clf;
set(gcf,'Position',framepos);
set(gcf,'Color','w');

%%% Standard deviation of all runs
for n=1:T
    s(n) = std(R(:,n));
    c1(n) = mean(R(:,n),1)+s(n);
    c2(n) = mean(R(:,n),1)-s(n);
end

%%% Plot Mean
pm = plot([1:T],mean(R(:,1:T),1),'color','k','linewidth',3);hold on;

%%% Plot shaded area of standard deviation
p = plot([1:T],c1); p.Color = 'none'; hold on; %xlim([1 T]);
p = plot([1:T],c2); p.Color = 'none';
T2 = [[1:T],fliplr([1:T])];
inBetween = [c1,fliplr(c2)];
pf = fill(T2,inBetween,[198/255,198/255,208/255]);
pf.EdgeColor = 'none';
pf.FaceAlpha = 0.5;

legend([pm pf],'mean','SD');


%%% Plot Reference case
p = plot([1:T],R(1,1:T),'color','k','linewidth',1); p.Color(4) = 0.5; hold on;


%%% Plot Hpyc (tan shades)
p = plot([1:T],R(5,1:T),'color',[236/255,221/255,154/255],'linewidth',1); p.Color(4) = 0.3; hold on;
p = plot([1:T],R(4,1:T),'color',[250/255,226/255,156/255],'linewidth',1); p.Color(4) = 0.3; 
p = plot([1:T],R(3,1:T),'color',[214/255,183/255,90/255],'linewidth',1); p.Color(4) = 0.3; 
p = plot([1:T],R(2,1:T),'color',[253/255,233/255,146/255],'linewidth',1); p.Color(4) = 0.3; 

%%% Plot H0 (pink shades)
p = plot([1:T],R(6,1:T),'color',[242/255,184/255,198/255],'linewidth',1); p.Color(4) = 0.3; hold on;
p = plot([1:T],R(7,1:T),'color',[246/255,154/255,205/255],'linewidth',1); p.Color(4) = 0.3; 
p = plot([1:T],R(8,1:T),'color',[253/255,92/255,168/255],'linewidth',1); p.Color(4) = 0.3; 
p = plot([1:T],R(9,1:T),'color',[225/255,21/255,132/255],'linewidth',1); p.Color(4) = 0.3; 

% Plot dS (red shades)
p = plot([1:T],R(10,1:T),'color',[1 0 0],'linewidth',1); p.Color(4) = 0.3;  hold on;
p = plot([1:T],R(11,1:T),'color',[0.8 0 0],'linewidth',1); p.Color(4) = 0.3; 
p = plot([1:T],R(12,1:T),'color',[0.6 0 0],'linewidth',1); p.Color(4) = 0.3; 
p = plot([1:T],R(13,1:T),'color',[0.4 0 0],'linewidth',1); p.Color(4) = 0.3; 

% Plot Tmax (purple shades)
p = plot([1:T],R(14,1:T),'color',[189/255,147/255,211/255],'linewidth',1); p.Color(4) = 0.7;  hold on;
p = plot([1:T],R(15,1:T),'color',[164/255,94/255,229/255],'linewidth',1); p.Color(4) = 0.3; 
p = plot([1:T],R(16,1:T),'color',[122/255,74/255,136/255],'linewidth',1); p.Color(4) = 0.3; 
p = plot([1:T],R(17,1:T),'color',[102/255,48/255,71/255],'linewidth',1); p.Color(4) = 0.7; 

% Plot Tatm (blue shades)
p = plot([1:T],R(18,1:T),'color',[0 0 1],'linewidth',2); p.Color(4) = 0.5;  hold on;
p = plot([1:T],R(19,1:T),'color',[0 0 0.73],'linewidth',1); p.Color(4) = 0.3; 
p = plot([1:T],R(20,1:T),'color',[0 0 0.47],'linewidth',1); p.Color(4) = 0.3; 
p = plot([1:T],R(21,1:T),'color',[0 0 0.3],'linewidth',2); p.Color(4) = 0.8; 

% Plot W (green shades)
p = plot([1:T],R(22,1:T),'color',[0.50,0.95,0.70],'linewidth',2); p.Color(4) = 0.5;  hold on;
p = plot([1:T],R(23,1:T),'color',[3/255,192/255,74/255],'linewidth',1); p.Color(4) = 0.3; 
p = plot([1:T],R(24,1:T),'color',[2/255,137/255,16/255],'linewidth',1); p.Color(4) = 0.3; 
p = plot([1:T],R(25,1:T),'color',[35/255,79/255,30/255],'linewidth',2); p.Color(4) = 0.7; 

% Plot rb (orange shades)
p = plot([1:T],R(26,1:T),'color',[253/255,174/255,29/255],'linewidth',1); p.Color(4) = 0.3;  hold on;
p = plot([1:T],R(27,1:T),'color',[237/255,112/255,20/255],'linewidth',1); p.Color(4) = 0.3; 
p = plot([1:T],R(28,1:T),'color',[201/255,92/255,10/255],'linewidth',1); p.Color(4) = 0.3; 
p = plot([1:T],R(29,1:T),'color',[137/255,49/255,1/255],'linewidth',1); p.Color(4) = 0.3; 


ax = gca;
ax.Box = 'on';
ax.YGrid='on';
ax.XGrid='on';
ax.FontName = 'Times New Roman';
ax.FontSize = 24;
xlim([1 30])