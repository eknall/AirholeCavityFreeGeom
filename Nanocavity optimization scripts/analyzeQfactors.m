clear all; close all; clc
%% Assemble loss anaylsis results
datLoc = 'D:\Files\FDTD simulations\Optimization runs\new design trial\';
files = dir([datLoc,'*.mat']);
data = [];
clear ds
for i = 1:length(files)
    load([datLoc,files(i).name]);
    data(i,:)= [ds.P.nholes ds.P.ndef ds.Qt ds.Qsc ds.Qwvg ds.Trans];
    clear ds
end

tmp = zeros(size(data));
[tmp(:,1) srtind] = sort(data(:,1));
for j = 1:length(data(1,:))
    if j ~= 1
        tmp(:,j) = data(srtind,j);
    end
end

%% Plot partial Q factors
figure;
set(gcf,'position',[969 249 944 1068-300])
axes('position',[0.1 0.1 0.85 0.85])
hold on;
box on

for k = 7
    pltind = find(tmp(:,2) == k);
    plot(tmp(pltind,1),(tmp(pltind,4)),'-.b','linewidth',1.5);
    plot(tmp(pltind,1),(tmp(pltind,5)),'-.r','linewidth',1.5);
    plot(tmp(pltind,1),(tmp(pltind,3)),'-k','linewidth',2.5);
end
ylim([1e2 1e7]);
xlim([7 17]);%max(tmp(pltind,1))
set(gca,'yscale','log','layer','top')

% title('n_{def} = 7','fontsize',16)
legend('Q_{sc}','Q_{wvg}','Q_{total}','location','northwest')
set(gca,'fontsize',14)
xlabel('# of mirror segments','fontsize',16)
ylabel('Q-factor','color','k','fontsize',16)

% Obtain the tick mark locations
xtick = 1:20;
% Obtain the limits of the y axis
ylimit = get(gca,'Ylim');
% Create line data
xX = repmat(xtick,2,1);
xY = repmat(ylimit',1,size(xtick,2));
% Plot line data
plot(xX,xY,'linestyle',':','color',0.75*ones(1,3),'linewidth',0.25)

% Obtain the tick mark locations
ytick = get(gca,'YTick');
yticks = [];
for i = 1:length(ytick)-1
    yticks = [yticks,linspace(ytick(i),ytick(i+1),10)];
end
        
% Obtain the limits of the y axis
xlimit = [1 20];
% Create line data
yY = repmat(yticks,2,1);
yX = repmat(xlimit',1,size(yticks,2));
% Plot line data
plot(yX,yY,'linestyle',':','color',0.75*ones(1,3),'linewidth',0.25)
% 
% annotation('textarrow',[0.825 0.825],[0.9 0.625], ...
%     'string','n_{def}','fontsize',16,'linewidth',2,'headstyle','cback3')
hold off


axes('position',[0.55 0.15 0.375 0.4]);
hold on;
box on

for k = 7
    pltind = find(tmp(:,2) == k);
    plot(tmp(pltind,1),(tmp(pltind,4))/10,'-.ob','linewidth',1.,'markersize',8);
    plot(tmp(pltind,1),(tmp(pltind,3)),'-ok','linewidth',2.,'markersize',8);
end

set(gca,'yscale','log','layer','top','xtick',7:9)
ylim([1.8e3 6.2e4]);
xlim([6.8 10.2]);
set(gca,'fontsize',14)
legend('Q_{sc} (/10)','Q_{total}','location','northwest')

% Plot line data
plot(xX,xY,'linestyle',':','color',0.75*ones(1,3),'linewidth',0.25)
plot(yX,yY,'linestyle',':','color',0.75*ones(1,3),'linewidth',0.25)



