% Draw 1/8 of the nanobeam cavity (using all available symmetry planes).
% The structure P is assumed to have the following fields, which define the
% geometry of the nanobeam cavity:
%
% P.a = lattice constant;
% P.wid = beam width; P.th = beam thickness;
% P.hhL = hole height Left side;
% P.hwL = hole width Left side;
% P.hhR = hole height Right side;
% P.hwR = hole width Right side;
% P.nholes = # of holes in 1/2 beam length;
% P.ndef = # of holes in 1/2 the defect region;
% P.maxdef = max. fractional increase in spacing/hole size in the defect
% P.oblong = hw is scaled by alpha^(1+oblong), hh by alpha^(1-oblong), e.g.
%   oblong == 1 results in constant height but changing width in defect;
%
% M.J. Burek, 08/16

function plotfieldprofilesXY(P,F,indx3d,QVdat)
%% Geometry parameters
hh = P.geom(:,1);
hw = P.geom(:,2);
xpos = P.geom(:,3);
ypos = P.geom(:,4);

%% plot3 field profile (w/ geometry) in the xy-plane
figure; set(gcf,'position',[9 1108-500 913 500])

ax = axes('position',[0.1 0.1 0.85 0.350]); hold(ax,'on')
surf(ax,F.xv,F.yv,real(F.Eden_z)'/max(max(real(F.Eden_z))));
text(-P.beamLen/2.1,1.25*P.w,max(max(real(F.Eden_z))),'|E|^2','color','w','fontname','times new roman', ...
    'fontangle','italic','fontsize',18)

bx = axes('position',[0.1 0.5 0.85 0.350]); hold(bx,'on')
surf(bx,F.xv,F.yv,F.Ey_z'/max(max(abs(F.Ey_z))));

text(-P.beamLen/2.1,1.25*P.w,max(max(abs(F.Ey_z))),'E_y','color','k','fontname','times new roman', ...
    'fontangle','italic','fontsize',18)

hdl = ax;
color = 'w';
for i = 1:2
    if i == 2
        hdl = bx;
        color = 'k';
    end
    
    daspect(hdl,[1,1,1])
    grid(hdl,'off')
    shading(hdl,'interp')
    view(2)
    
    %waveguide overlay
    
    %waveguide overlay
    indx.in = indx3d.index_x;
    len = length(indx.in(1,1,:));
    x=indx3d.x;
    y = indx3d.y;
    
    %Plot top of structure
    %find top
    l = len;
    
    while(mean(mean(indx.in(:,:,l))) <= 1)
        l = l - 1;
    end
    indz = squeeze(real(indx.in(:,:,l)))';
    contour3(hdl,x,y,indz,2.6,color);
    
%     
%     plot3(hdl,x,y,indz,color,'linewidth',1)
%     plot3(hdl,-x,y,indz,color,'linewidth',1)
    
%     plot3(hdl,1/2*[-P.beamLen P.beamLen],1/2*[-P.w -P.w],ones(2),color,'linewidth',1)
%     plot3(hdl,1/2*[-P.beamLen P.beamLen],1/2*[P.w P.w],ones(2),color,'linewidth',1)
    
    %airhole overlay
%     for j=1:length(hw);
%         xdat = linspace(-hh(j)/2,hh(j)/2,100)';
%         ydatu(:,j) = sqrt(1-linspace(-hh(j)/2,hh(j)/2,100).^2/(hh(j)/2)^2)*hw(j)/2 + ypos(j);
%         ydatd(:,j) = -sqrt(1-linspace(-hh(j)/2,hh(j)/2,100).^2/(hh(j)/2)^2)*hw(j)/2 + ypos(j);
%                 
%         plot3(hdl,[xdat + xpos(j),xdat + xpos(j)], ...
%             [ydatu(:,j),ydatd(:,j)], ...
%             ones(size(xdat)),'color',color,'linewidth',1);        
%     end
    plot3(hdl,[0,0],2*P.w*[-1,1], ...
        [1,1],'r','linestyle',':','linewidth',0.5);
    
    set(hdl,'xtick',[],'ytick',[],'xticklabel','','yticklabel','');
    hold(hdl,'off');
    box(hdl,'off');
    
    daspect(hdl,[1,1,1])
    
    xlim(hdl,[-P.beamLen/2 P.beamLen/2]);
    hlim = axis(hdl);
    ylim(hdl,[hlim(1)/5 hlim(2)/5]);
end

linkaxes([ax,bx],'x')

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

title(ax, ...
    {['\lambda_o= ',num2str(QVdat.lambda*1e9,'%.1f'),'nm ', ...
    'V_{mode}= ',num2str(QVdat.Vmode,'%.2f'), '(\lambda/n)^3 ', ...
    'Q_t= ',num2str(roundn(QVdat.Qt,3),'%.0f'), ' ',]; ...
    ['Q_{sc}= ',num2str(roundn(QVdat.Qsc,3),'%.0f'), ' ', ...
    'Q_{wvg}= ',num2str(roundn(QVdat.Qwvg,3),'%.0f'), ' ', ...
    'T = ',num2str(100*QVdat.Trans,'%.0f'),'%']}, ...
    'fontname','arial','fontsize',12);





