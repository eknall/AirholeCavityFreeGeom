clear all; close all
%% example plot

omega0 = 1;
omega = (0.99:1e-8:1.01)*omega0;

Qout = 1e6;
Qin = fliplr(logspace(4,6,7));
Qrad = 1e5;

cmap = jet(length(Qin))./2;
linstyl = {'-','--','-.',':','-.','--','-'};

tauOut = 2*Qout/omega0;
tauRad = 2*Qrad/omega0;

T = @(omega,tauIn) (4/(tauIn*tauOut))./((omega-omega0).^2 + (1/tauIn + 1/tauOut + 1/tauRad)^2);
R = @(omega,tauIn) ((omega-omega0).^2 + (1/tauRad - 1/tauIn + 1/tauOut)^2)./((omega-omega0).^2 + (1/tauIn + 1/tauOut + 1/tauRad)^2);

figure;
set(gcf,'position',[10,49,943*0.75,1068]);
ax=axes('position',[0.125 0.325 0.8 0.6]);
hold(ax,'on')

for k = 1:length(Qin);
    tauIn = 2*Qin(k)/omega0;
    
    set(ax,'linestyle',linstyl(k));  
    plot(ax,omega/omega0,R(omega,tauIn),'linewidth',2.5,'Color',cmap(k,:))
    text(1.0000125+5e-6*(k),0.5-0.05*(k),['Q_{in} = ',num2str(Qin(k),'%.1e')],'fontsize',18,'Color',cmap(k,:))
end
title(ax,['Q_{rad} = ',num2str(Qrad,'%.1e'),', Q_{out} = ',num2str(Qout,'%.1e')],'fontsize',20)
ylabel(ax,'R(\omega)','fontsize',20)
xlim(ax,[0.9999 1.0001])
ylim(ax,[0 1])
set(ax,'fontsize',18,'ytick',0:0.2:1,'xticklabel',{})
box on
hold(ax,'off')


bx=axes('position',[0.125 0.1 0.8 0.2]);
hold(bx,'on')
for k = 1:length(Qin);
    tauIn = 2*Qin(k)/omega0;    
    set(bx,'linestyle',linstyl(k));  
    plot(bx,omega/omega0,T(omega,tauIn),'k','linewidth',2.5,'Color',cmap(k,:))
end

xlim(bx,[0.9999 1.0001])
ylim(bx,[0 0.21])

set(bx,'fontsize',18,'xtick',0.999:0.001:1.001,'ytick',[0,0.1,0.2])
xlabel(bx,'\omega/\omega_o','fontsize',20)
ylabel(bx,'T(\omega)','fontsize',20)
box on


%% single plot

figure;
set(gcf,'position',[10,49,943*0.75,1068]);
ax=axes('position',[0.125 0.325 0.8 0.55]);
hold(ax,'on')

omega0 = 2*pi*196e12;
omega = (0.99:1e-8:1.01)*omega0;

Qout = 2.3e5;
Qin = 2.3e5;
Qrad = 1.8e5;
Qtot = (1/Qrad + 1/Qin + 1/Qout)^-1;

tauOut = 2*Qout/(omega0/2/pi);
tauRad = 2*Qrad/(omega0/2/pi);
tauIn = 2*Qin/(omega0/2/pi);

T = @(omega,tauIn) (4/(tauIn*tauOut))./((omega-omega0).^2 + (1/tauIn + 1/tauOut + 1/tauRad)^2);
R = @(omega,tauIn) ((omega-omega0).^2 + (1/tauRad - 1/tauIn + 1/tauOut)^2)./((omega-omega0).^2 + (1/tauIn + 1/tauOut + 1/tauRad)^2);

kappaRad = 1e-9/tauRad;
kappaIn = 1e-9/tauIn;
kappaOut = 1e-9/tauOut;
kappaTot = kappaIn + kappaRad + kappaOut;

plot(ax,omega/omega0,R(omega,tauIn),'k','linewidth',2.5)
text(1.00000125,0.5,['Q_{tot} = ',num2str(Qtot,'%.1e')],'fontsize',20)
text(1.00000125,0.5-0.075,['\kappa_{tot} = ',num2str(kappaTot,'%.2f'),'GHz'],'fontsize',20)

    
title(ax,{['Q_{rad} = ',num2str(Qrad,'%.1e'),', Q_{in} = ',num2str(Qin,'%.1e'),', Q_{out} = ',num2str(Qout,'%.1e')] ...
    ; ['\kappa_{rad} = ',num2str(kappaRad,'%.2f'),'GHz',', \kappa_{in} = ',num2str(kappaIn,'%.2f'),'GHz', ...
    ', \kappa_{out} = ',num2str(kappaOut,'%.2f'),'GHz']},'fontsize',20)
ylabel(ax,'R(\omega)','fontsize',20)
xlim(ax,[0.99999 1.00001])
ylim(ax,[0 1])
set(ax,'fontsize',18,'ytick',0:0.2:1,'xticklabel',{})
box on
hold(ax,'off')


bx=axes('position',[0.125 0.1 0.8 0.2]);
hold(bx,'on')
tauIn = 2*Qin/omega0;    
plot(bx,omega/omega0,T(omega,tauIn),'k','linewidth',2.5)

xlim(bx,[0.99999 1.00001])
ylim(bx,[0 1.01])

set(bx,'fontsize',18,'xtick',0.999:0.001:1.001,'ytick',0:0.2:1)
xlabel(bx,'\omega/\omega_o','fontsize',20)
ylabel(bx,'T(\omega)','fontsize',20)
box on
