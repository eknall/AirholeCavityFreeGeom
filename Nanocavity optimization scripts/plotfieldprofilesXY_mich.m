
% The structure P is assumed to have the following fields, which define the
% geometry of the nanobeam cavity:
%
% P.a = lattice constant;
% P.wid = beam width; 
% P.theta = beam half angle;


%
% M.J. Burek, 08/16
% Modified by M. Chalupnik 2/6/19

function plotfieldprofilesXY_mich(P,F,indx3d,QVdat)
%% Geometry parameters

bcutoff = 0.6;

%% plot3 field profile (w/ geometry) in the xy-plane
figure; set(gcf,'position',[9 1108-500 913 600])

%ax is E^2 (intensity)
ax = axes('position',[0.1 0.1 0.85 0.350]); 
hold(ax,'on')

xlim(ax,[-P.beamLen*bcutoff/2 P.beamLen*bcutoff/2]);
bx = axes('position',[0.1 0.55 0.85 0.350]); hold(bx,'on')


hdl = ax;
color = 'w';
colormap(ax, jet);
caxis(ax, [0,1]);
colormap(bx,bwr);
caxis(bx, [-1,1]);

for i = 1:2
    if i == 2
        hdl = bx;
        color = 'k';
    end
    

    
       
   
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
    contour3(hdl, -x, y, indz, 2.6, color);
    
    if i ==1
        Eden =  F.n_x.^2 .* abs(F.Ex).^2 + F.n_y.^2 .* abs(F.Ey).^2 + F.n_z.^2 .* abs(F.Ez).^2;
        
        Edenx = squeeze(Eden(:, :,l))/max(max(max(Eden)));
        surf(ax,F.xv,F.yv,Edenx');
            daspect(hdl,[1,1,1])
    grid(hdl,'off')
    shading(hdl,'interp')
    view(2)
%         text(ax,-0.9*P.w,0.65*P.w,max(max(Edenx)),'|E|^2','color','w','fontname','times new roman', ...
%             'fontangle','italic','fontsize',18)
    else
        Ey = squeeze(real(F.Ey( :, :, l)));
        
        Ey_x = Ey / max(max(max(abs(Ey))));
        surf(bx,F.xv,F.yv, Ey_x');
            daspect(hdl,[1,1,1])
    grid(hdl,'off')
    shading(hdl,'interp')
    view(2)
%         text(bx,-0.9*P.w,0.65*P.w,max(max(Ey_x)),'E_y','color','k','fontname','times new roman', ...
%             'fontangle','italic','fontsize',18)
    end
   

    set(hdl,'xtick',[],'ytick',[],'xticklabel','','yticklabel','');
   % hold(hdl,'off');
    box(hdl,'off');
    
    daspect(hdl,[1,1,1])
    
    xlim(hdl,[-P.beamLen*bcutoff/2 P.beamLen*bcutoff/2]);
    hlim = axis(hdl);
    ylim(hdl,[hlim(1)/10 hlim(2)/10]); %how why this in particular?
    
    hold(hdl,'off');
end

linkaxes([ax,bx],'x') %not sure about this

title(bx, ...
    {['\theta = ',num2str(P.theta),'^o; ', 'w = ',num2str(P.w*1e9,'%.0f'),'nm; ']; ...
    ['Mirror Geom: a = ',num2str(P.a(1)*1e9,'%.0f'),'nm; ',... 
    '\alpha= ',num2str(P.alpha(1),'%.0f'),'^o; ', ...
    '\phi= ',num2str(P.phi(1),'%.0f'),'^o; ', ...
    'orad= ',num2str(P.orad(1)*1e9,'%.0f'),'nm; ', ...
    'irad= ',num2str(P.irad(1)*1e9,'%.0f'),'nm; ']; ...
        ['Cavity Geom: a = ',num2str(P.a(2)*1e9,'%.0f'),'nm; ',... 
    '\alpha= ',num2str(P.alpha(2),'%.0f'),'^o; ', ...
    '\phi= ',num2str(P.phi(2),'%.0f'),'^o; ', ...
    'orad= ',num2str(P.orad(2)*1e9,'%.0f'),'nm; ', ...
    'irad= ',num2str(P.irad(2)*1e9,'%.0f'),'nm; ']; ...
    ['n_{cell}= ',num2str(P.ncell),'; ', ...
    'n_{half mir cell}= ',num2str(P.nmircellh),' ']}, ...
    'fontname','arial','fontsize',14)

title(ax, ...
    {['\nu_o =', num2str((3*10^8 /QVdat.lambda) *10^-12,'%.1f'), 'THz, \lambda_o= ',num2str(QVdat.lambda*1e9,'%.1f'),'nm ', ...
    'V_{mode}= ',num2str(QVdat.Vmode,'%.2f'), '(\lambda/n)^3 ', ...
    'Q_t= ',num2str(roundn(QVdat.Qt,1),'%.0f'), ' ',]; ...
    ['Q_{sc}= ',num2str(roundn(QVdat.Qsc,1),'%.0f'), ' ', ...
    'Q_{wvg}= ',num2str(roundn(QVdat.Qwvg,1),'%.0f'), ' ', ...
    'T = ',num2str(100*QVdat.Trans,'%.0f'),'%']}, ...
    'fontname','arial','fontsize',14);





