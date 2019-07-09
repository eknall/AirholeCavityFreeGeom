% Draw 1/8 of the nanobeam cavity (using all available symmetry planes).
% The structure P is assumed to have the following fields, which define the
% geometry of the nanobeam cavity:
%
% P.a = lattice constant;
% P.wid = beam width; P.th = beam thickness;
% P.hh = hole height;
% P.hw = hole width;
% P.nholes = # of holes in 1/2 beam length;
% P.ndef = # of holes in 1/2 the defect region;
% P.maxdef = max. fractional increase in spacing/hole size in the defect
% P.oblong = hw is scaled by alpha^(1+oblong), hh by alpha^(1-oblong), e.g.
%   oblong == 1 results in constant height but changing width in defect;
%
% M.J. Burek, 08/16

 function plotfieldprofilesYZ(P,F,QVdat)
%% Geometry parameters
asym = P.asym;

%% plot3 field profile (w/ geometry) in the xy-plane
figure; set(gcf,'position',[9 1108-500 913 500])

ax = axes('position',[0.1 0.1 0.85 0.250]); hold(ax,'on')
surf(ax,F.yv,F.zv,real(F.Eden_x)'/max(max(real(F.Eden_x))));
text(-0.9*P.w,0.65*P.w,min(min(real(F.Eden_z))),'|E|^2','color','w','fontname','times new roman', ...
    'fontangle','italic','fontsize',18)
text(-0.9*P.w,0.65*P.w,min(min(real(F.Eden_z))),'|E|^2','color','w','fontname','times new roman', ...
    'fontangle','italic','fontsize',18)

bx = axes('position',[0.1 0.5 0.85 0.250]); hold(bx,'on')
surf(bx,F.yv,F.zv,F.Ey_x'/max(max(abs(F.Ey_x))));

text(-0.9*P.w,0.65*P.w,max(max(F.Ey_x)),'E_y','color','w','fontname','times new roman', ...
    'fontangle','italic','fontsize',18)
text(-0.9*P.w,0.65*P.w,min(min(F.Ey_x)),'E_y','color','w','fontname','times new roman', ...
    'fontangle','italic','fontsize',18)

hdl = ax;
for i = 1:2
    if i == 2
        hdl = bx;
    end
    
    daspect(hdl,[1,1,1])
    grid(hdl,'off')
    shading(hdl,'interp')
    view(2)
    
    %waveguide cross-section overlay
    plot3(hdl,1/2*[-P.w P.w], 1/2*[P.th P.th], ones(2),'w','linewidth',1)
    plot3(hdl,1/2*[-P.w 2*asym], 1/2*[P.th, -P.th], ones(2),'w','linewidth',1)
    plot3(hdl,1/2*[2*asym P.w], 1/2*[-P.th, P.th], ones(2),'w','linewidth',1)
    
    plot3(hdl,[0,0],2*P.w*[-1,1], ...
        [1,1],'r','linestyle',':','linewidth',0.5);
    
    if asym > 0
        plot3(hdl,[asym,asym],2*P.w*[-1,1], ...
            [1,1],'w','linestyle',':','linewidth',0.5);
    end
    
    set(hdl,'xtick',[],'ytick',[],'xticklabel','','yticklabel','');
    hold(hdl,'off');
    box(hdl,'off');
    
    daspect(hdl,[1,1,1])

    xlim(hdl,[-P.w P.w]);
    hlim = axis(hdl);
    ylim(hdl,[hlim(1) hlim(2)]);    
end

linkaxes([ax,bx],'x')

title(bx, ...
    {['\theta = ',num2str(P.theta),'^o ',...
    'aL = ',num2str(P.aL*1e9,'%.0f'),'nm ',...
    'aR = ',num2str(P.aR*1e9,'%.0f'),'nm ',...
    'w = ',num2str(P.w*1e9,'%.0f'),'nm '];...
    ['h_xL= ',num2str(P.hhL*1e9,'%.0f'),'nm ', ...
    'h_yL= ',num2str(P.hwL*1e9,'%.0f'),'nm ', ...
    'h_xR= ',num2str(P.hhR*1e9,'%.0f'),'nm ', ...
    'h_yR= ',num2str(P.hwR*1e9,'%.0f'),'nm ']; ...
    ['n_{hole}= ',num2str(P.nholes),' ', ...
    'n_{def}= ',num2str(P.ndef),' ', ...
    'd_{max}= ',num2str(P.maxdef),' ', ...
    'oblong = ',num2str(P.oblong)]}, ...
    'fontname','arial','fontsize',14)

title(ax, ...
    {['\lambda_o= ',num2str(QVdat.lambda*1e9,'%.1f'),'nm ', ...
    'V_{mode}= ',num2str(QVdat.Vmode,'%.2f'), '(\lambda/n)^3 ', ...
    'Q_t= ',num2str(roundn(QVdat.Qt,3),'%.0f'), ' ',]; ...
    ['Q_{sc}= ',num2str(roundn(QVdat.Qsc,3),'%.0f'), ' ', ...
    'Q_{wvg}= ',num2str(roundn(QVdat.Qwvg,3),'%.0f'), ' ', ...
    'T = ',num2str(100*QVdat.Trans,'%.0f'),'%']}, ...
    'fontname','arial','fontsize',14);





