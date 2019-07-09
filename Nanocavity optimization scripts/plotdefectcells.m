% Draw the nanobeam cavity (using all available symmetry planes).
% The structure P is assumed to have the following fields, which define the
% geometry of the nanobeam cavity:
%
% P.aL = nominal lattice constant Left Side
% P.aR = nominal lattice constant Right Side
% P.w = beam width; 
% P.th = beam thickness;
% P.hhL = hole height Left Side;
% P.hwL = hole width Left Side;;
% P.hhR = hole height Right Side;
% P.hwR = hole width Right Side;
% P.nholes = # of holes in 1/2 beam length;
% P.ndef = # of holes in 1/2 the defect region;
% P.maxdef = max. fractional increase in spacing/hole size in the defect
% P.oblong = hw is scaled by alpha^(1+oblong), hh by alpha^(1-oblong), e.g.
%   oblong == 1 results in constant height but changing width in defect;
%
% M.J. Burek, 08/16
% E.N. Knall, 04/19

function plotdefectcells(P)
%% Geometry parameters

a = P.ahole;

hh = P.geom(:,1);
hw = P.geom(:,2);
xpos = P.geom(:,3);
ypos = P.geom(:,4);
asym = P.asym;

%% Plot defect region
figure; set(gcf,'position',[9 1108-500 913 500])
ax = axes('position',[0.1 0.1 0.85 0.4]);
hold(ax,'on')

% assemble center positions for each mirror segment
xdat = zeros(1,(length(xpos)-1));
a_hole = zeros(1,(length(a)-1));
a_hole(1:floor(end/2)) = a(1:end/2-1);
a_hole(floor(end/2)+1:end) = a(end/2+1:end);
%xdat defines lattice constant marker positions
for j=1:length(xdat)
    xdat(j) = (xpos(j+1)+xpos(j))/2;
end
%it's necessary to flip the position of the holes because
%we flip the order of the holes in runNanobeamCavity
%xdat = fliplr(xdat);
% for j=2:length(xdat)
%     xdat(j) = -xdat(j)+(a(j-1)+a(j))/2;
% end
hh_hole = hh';
hw_hole = hw';

hdl = plot(ax,xdat,ones(length(xdat)),'--r','LineWidth',1.5);
%arbitrarily choosing to normalize hole heights by aR 
hdl(2) = plot(ax,xdat,a_hole/P.aR,'-ks','LineWidth',1, ...
    'MarkerSize',6,'Markerfacecolor','k');

% assemble center positions for each hole position
hdl(3) = plot(ax,xpos,hh_hole/P.aR,'-s','LineWidth',1, ...
    'color',[0 102 51]/255,'MarkerSize',6,'Markerfacecolor',[0 102 51]/255);
hdl(4) = plot(ax,xpos,hw_hole/P.aR,'-bs','LineWidth',1, ...
    'MarkerSize',6,'Markerfacecolor','b');

axis(ax,'tight')
axl = axis(ax);
ylim(ax,[0, 0.05+axl(4)]);
xlabel(ax,'mirror segment #','FontName','arial','FontSize',14);
ylabel(ax,'a_{nominal}','FontName','arial','FontSize',16);

legend(hdl(2:end),'a','h_x','h_y','location','southeast');
set(ax,'FontName','arial','FontSize',14,'xtick',[],'xticklabel','');

hold(ax,'off');
box(ax,'on');


%% Plot geometry in the xy-plane

bx = axes('position',[0.1 0.475 0.85 0.4]);
hold(bx,'on')

%waveguide overlay
plot(bx,1/2*[-P.beamLen P.beamLen],1/2*[-P.w -P.w],'k','linewidth',1.)
plot(bx,1/2*[-P.beamLen P.beamLen],1/2*[P.w P.w],'k','linewidth',1)

%airhole overlay
if isfield(P, 'wvgmir')   
    if isvector(P.wvgmir)
        wvgmirlength = P.wvgmir;
    else
        wvgmirlength = length(P.wvgmir);
    end    
else
    wvgmirlength = 0;
end

for j=1:length(hw)
    xdat = linspace(-hh(j)/2,hh(j)/2,100)';
    ydatu(:,j) = sqrt(1-linspace(-hh(j)/2,hh(j)/2,100).^2/(hh(j)/2)^2)*hw(j)/2 + ypos(j);
    ydatd(:,j) = -sqrt(1-linspace(-hh(j)/2,hh(j)/2,100).^2/(hh(j)/2)^2)*hw(j)/2 + ypos(j);
    
    
    plot(bx,[xdat + xpos(j),xdat + xpos(j)], ...
        [ydatu(:,j),ydatd(:,j)],'k','linewidth',1);
    
    
    if j > (length(hw)/2 - P.ndef) && j < (length(hw)/2 + P.ndef+1)
        plot(bx,[xpos(j),xpos(j)],2*P.w*[-1,1], ...
            'k','linestyle',':','linewidth',.5);
    elseif j == (length(hw)/2 - P.ndef) || j == (length(hw)/2 + P.ndef+1)
        plot(bx,[xpos(j),xpos(j)],2*P.w*[-1,1], ...
            'b','linestyle','--','linewidth',.5);
    elseif j == wvgmirlength+1 || j == (length(hw) - wvgmirlength)
        plot(bx,[xpos(j),xpos(j)],2*P.w*[-1,1], ...
            'b','linestyle','--','linewidth',.5);
    end
end
plot(bx,[0,0],2*P.w*[-1,1], ...
    'r','linestyle','--','linewidth',.5);
set(bx,'xtick',[],'ytick',[]);
hold(bx,'off');
box(bx,'on')
set(bx,'xcolor','w','ycolor','w');
daspect([1 1 1])
xlim(bx,[-P.beamLen/2 P.beamLen/2]);
bxl = axis(bx);
ylim(bx,[bxl(1)/5 bxl(2)/5]);
linkaxes([bx,ax],'x')

title(bx, ...
    {['\theta = ',num2str(P.theta),'^o ',...
    'aL = ',num2str(P.aL*1e9,'%.0f'),'nm ',...
    'aR = ',num2str(P.aR*1e9,'%.0f'),'nm ',...
    'w = ',num2str(P.w*1e9,'%.0f'),'nm ', ...
    'h_xL= ',num2str(P.hhL*1e9,'%.0f'),'nm ', ...
    'h_yL= ',num2str(P.hwL*1e9,'%.0f'),'nm '];...
    ['h_xR= ',num2str(P.hhR*1e9,'%.0f'),'nm ', ...
    'h_yR= ',num2str(P.hwR*1e9,'%.0f'),'nm ', ...
    'n_{hole}= ',num2str(P.nholes),' ', ...
    'n_{def}= ',num2str(P.ndef),' ', ...
    'd_{max}= ',num2str(P.maxdef),' ', ...
    'oblong = ',num2str(P.oblong)]}, ...
    'fontname','arial','fontsize',12)


%% Plot geometry in the yz-plane
if asym > 0
    cx = axes('position',[0.825 0.8 0.125 0.125]);
    hold(cx,'on')
    
    %waveguide overlay
    plot(cx,1/2*[-P.w P.w], 1/2*[P.th P.th],'k','linewidth',1)
    plot(cx,1/2*[-P.w -2*asym], 1/2*[P.th, -P.th],'k','linewidth',1)
    plot(cx,1/2*[-2*asym P.w], 1/2*[-P.th, P.th],'k','linewidth',1)
        
    plot(cx,[0,0],1.05*P.th*[-1,1], ...
        'r','linestyle','--','linewidth',.5);
    plot(cx,[-asym,-asym],1.05*P.th*[-1,1], ...
        'k','linestyle','--','linewidth',.5);
   
    set(cx,'xtick',[],'ytick',[]);
    hold(cx,'off');
    box(cx,'on')
    set(cx,'xcolor','w','ycolor','w');
    daspect([1 1 1])
    xlim(cx,0.65*[-P.w P.w]);
    ylim(cx,[-P.th P.th]);  
    text(0,-(P.th/2),[' ~',num2str(asym*1e9,'%.0f'),'nm'],'fontsize',10)
end





