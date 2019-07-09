% Draw the nanobeam cavity (using all available symmetry planes).
% The structure P is assumed to have the following fields, which define the
% geometry of the nanobeam cavity:
%
% P.a = lattice constant;
% P.w = beam width; P.th = beam thickness;
% P.hh = hole height;
% P.hw = hole width;
% P.nholes = # of holes in 1/2 beam length;
% P.ndef = # of holes in 1/2 the defect region;
% P.maxdef = max. fractional increase in spacing/hole size in the defect
% P.oblong = hw is scaled by alpha^(1+oblong), hh by alpha^(1-oblong), e.g.
%   oblong == 1 results in constant height but changing width in defect;
%
% M.J. Burek, 08/16
% E.N. Knall, 04/19

function P = createNanobeamCavity(P)

a = P.a;
wid = P.w;
thi = P.th;
hhL = P.hhL;
hwL = P.hwL;
hhR = P.hhR;
hwR = P.hwR;
nholes = P.nholes;
ndef = P.ndef;
maxdef = P.maxdef;
oblong = P.oblong;

%% Check that the parameters are physically reasonable
% Commented while converting to asymmetric design EK
% if wid-hw < 100e-9 
%     error('Hole width too large relative to beam width');
% elseif hw < 125e-9 || hh < 125e-9
%     error('Hole diameter too small for lithography constraints');
% end
% 
% if P.consthole == 1  
%     if a*(1-maxdef)-hh < 50e-9
%         error('Hole height too large relative to lattice constant');
%     end
% else   
%     if a*(1-maxdef)-hh*(1-maxdef)^(1-oblong) < 50e-9
%         error('Hole height too large relative to lattice constant');
%     elseif hh*(1-maxdef)^(1-oblong) < 125e-9 || hw*(1-maxdef)^(1+oblong) < 125e-9
%         error('Defect hole diameter too small for lithography constraints');
%     end   
% end

%% Assemble geometry
totNumHoles = 2*nholes;
lenWvgMir = 0;
if isfield(P, 'wvgmir')    
    if isvector(P.wvgmir)
        lenWvgMir = P.wvgmir;
    else
        lenWvgMir = length(P.wvgmir(:,1));
    end    
end
totNumHoles = totNumHoles + lenWvgMir*2;

a_hole = zeros(1,totNumHoles);
hh_hole = zeros(1,totNumHoles);
hw_hole = zeros(1,totNumHoles);
xpos = zeros(1,totNumHoles);

holeIndex = lenWvgMir+nholes+1;
hhTransSlope = (hhL-hhR)/(2*P.ndef - 1);
hwTransSlope = (hwL-hwR)/(2*P.ndef - 1);
for k = 1:nholes
    hhRTrans = hhR-(k-P.ndef)*hhTransSlope;
    hwRTrans = hwR-(k-P.ndef)*hwTransSlope;
    if k == 1
        a_hole(holeIndex) = a*(1-maxdef);
        %offset so dielectric is centered
        xpos(holeIndex-1) = -a_hole(holeIndex)/2;
        if P.consthole == 1
            hh_hole(holeIndex) = hhRTrans;
            hw_hole(holeIndex) = hwRTrans;
        else
            hh_hole(holeIndex) = hhRTrans*(1-maxdef)^(1-oblong);
            hw_hole(holeIndex) = hwRTrans*(1-maxdef)^(1+oblong);
        end
    elseif k <= P.ndef
        a_hole(holeIndex) = a*(1-maxdef* ...
            (2*((k-1)/ndef)^3-3*((k-1)/ndef)^2+1));
        if P.consthole == 1
            hh_hole(holeIndex) = hhRTrans;
            hw_hole(holeIndex) = hwRTrans;
        else
            hh_hole(holeIndex) = hhRTrans*(1-P.maxdef* ...
                (2*((k-1)/ndef)^3-3*((k-1)/ndef)^2+1))^(1-oblong);
            hw_hole(holeIndex) = P.hwRTrans*(1-P.maxdef* ...
                (2*((k-1)/ndef)^3-3*((k-1)/ndef)^2+1))^(1+oblong);
        end
    else
        a_hole(holeIndex) = a;
        hh_hole(holeIndex) = hhR;
        hw_hole(holeIndex) = hwR;
    end
    xpos(holeIndex) = xpos(holeIndex-1)+a_hole(holeIndex);
    holeIndex = holeIndex + 1;
end

% append wvg-mirror taper if defined
if isfield(P, 'wvgmir')    
    if isvector(P.wvgmir)
        wvgmir = P.wvgmir+1; %eliminate the first mirror hole and null hole        
        for k = 2:wvgmir
            a_hole(holeIndex) = a*(1-((k-1)/wvgmir)^2);
            hh_hole(holeIndex) = hhR*(1-((k-1)/wvgmir)^2);
            hw_hole(holeIndex) = hwR*(1-((k-1)/wvgmir)^2);
            xpos(holeIndex) = xpos(holeIndex-1)+a_hole(holeIndex);
            holeIndex = holeIndex + 1;
        end        
    else
        for k = 1:length(P.wvgmir(:,1))
            a_hole(holeIndex) = P.wvgmir(k,1);
            hh_hole(holeIndex) = P.wvgmir(k,2);
            hw_hole(holeIndex) = P.wvgmir(k,3);
            xpos(holeIndex) = xpos(holeIndex-1)+a_hole(holeIndex)/2;
            holeIndex = holeIndex + 1;
        end
    end    
end

holeIndex = lenWvgMir+nholes;
for k = 1:nholes
    hhLTrans = hhL+(k-P.ndef)*hhTransSlope;
    hwLTrans = hwL+(k-P.ndef)*hwTransSlope;
    if k == 1
        a_hole(holeIndex) = a*(1-maxdef);
        if P.consthole == 1
            hh_hole(holeIndex) = hhLTrans;
            hw_hole(holeIndex) = hwLTrans;
        else
            hh_hole(holeIndex) = hhLTrans*(1-maxdef)^(1-oblong);
            hw_hole(holeIndex) = hwLTrans*(1-maxdef)^(1+oblong);
        end
    elseif k <= P.ndef
        a_hole(holeIndex) = a*(1-maxdef* ...
            (2*((k-1)/ndef)^3-3*((k-1)/ndef)^2+1));
        if P.consthole == 1
            hh_hole(holeIndex) = hhLTrans;
            hw_hole(holeIndex) = hwLTrans;
        else
            hh_hole(holeIndex) = hhLTrans*(1-P.maxdef* ...
                (2*((k-1)/ndef)^3-3*((k-1)/ndef)^2+1))^(1-oblong);
            hw_hole(holeIndex) = P.hwLTrans*(1-P.maxdef* ...
                (2*((k-1)/ndef)^3-3*((k-1)/ndef)^2+1))^(1+oblong);
        end
    else
        a_hole(holeIndex) = a;
        hh_hole(holeIndex) = hhL;
        hw_hole(holeIndex) = hwL;
    end
    xpos(holeIndex) = xpos(holeIndex+1)-a_hole(holeIndex);
    holeIndex = holeIndex - 1;
end

% append wvg-mirror taper if defined
if isfield(P, 'wvgmir')    
    if isvector(P.wvgmir)
        wvgmir = P.wvgmir+1; %eliminate the first mirror hole and null hole        
        for k = 2:wvgmir
            a_hole(holeIndex) = a*(1-((k-1)/wvgmir)^2);
            hh_hole(holeIndex) = hhL*(1-((k-1)/wvgmir)^2);
            hw_hole(holeIndex) = hwL*(1-((k-1)/wvgmir)^2);
            xpos(holeIndex) = xpos(holeIndex+1)-a_hole(holeIndex);
            holeIndex = holeIndex - 1;
        end        
    else
        for k = 1:length(P.wvgmir(:,1))
            a_hole(holeIndex) = P.wvgmir(k,1);
            hh_hole(holeIndex) = P.wvgmir(k,2);
            hw_hole(holeIndex) = P.wvgmir(k,3);
            xpos(holeIndex) = xpos(holeIndex+1)-a_hole(holeIndex);
            holeIndex = holeIndex - 1;
        end
    end    
end


% assemble air hole center positions
ypos = zeros(size(a_hole));

% output final geometry
P.ahole = a_hole';
P.beamLen = sum(a_hole)+a;
P.geom = [hh_hole' hw_hole' xpos' ypos'];

