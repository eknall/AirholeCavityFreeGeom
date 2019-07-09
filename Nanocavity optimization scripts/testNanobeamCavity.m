clear all; close all; clc

P.aL = 280e-9;                          % 260nominal lattice constant Left Side
P.aR = 250e-9;                          % 280nominal lattice constant Right Side
P.w = 495e-9;                           % beam width 
P.theta = 50;                           % etch angle in degrees
P.th = P.w/(2*tan(P.theta*pi/180));     % beam thickness

P.hhL = 150e-9;                          % nominal hole height Left Side
P.hwL = 150e-9;                          % nominal hole width Left Side
P.hhR = 160e-9;                          % nominal hole height Right Side
P.hwR = 160e-9;                          % nominal hole width Right Side

P.nholes = 12;                          % # holes in 1/2 beam length
P.ndef = 5;                             % # of holes in 1/2 defect region

% OPTIONAL
P.wvgmir = 5;                               % # of segments of appended linear wvg-mir taper
% P.wvgmir = 1e-9*[274.98	151.86	160.82  % custom array of end wvg-mir taper ([a;hh;hw])
%                  263.00	137.21	143.69
%                  248.70	119.72	123.24
%                  236.72	105.07	106.11
%                  231.70	98.93	98.93];

%P.maxdef = 0.098264;                        % defect percentage
P.maxdef = 0.10006;                    % defect percentage
P.consthole = 1;                        % 1 if hole size is held constant
P.oblong = 1.5*(1-P.consthole);           % oblong parameter (zero if holes are not changed)

P.nbeam = 2.4028;                       % index of refraction in material
P.lambda = 737e-9;                      % target optical wavelength 740

%Disorder
P.stdDev = 0.0*[P.hhL,P.hwL];             % standard deviation of hole dimensions (hh,hw)
P.stdDevPos = 0.0*P.aR;                  % standard deviation of hole positions
P.asym = 0.0*P.w;                       % cross-section asymmetry (target y-offset in bottom apex position)

P.storeFDTD = 1;                        % 1 to keep solved FDTD file


%% required FDTD file locations

% directory path for location of FDTD files
files.path = [pwd,'/Nanocavity optimization scripts/FDTD scripts/']; % Update for location of script files
files.template = [files.path,'FDTDtemplate.fsp']; % template file containing all preset FDTD elements/parameters
files.build = [files.path,'script_buildNanobeamCavity.lsf']; % Lumerical script to build geomtry
files.runsim = [files.path,'script_runNanobeamCavity.lsf']; % Lumerical script to run simulation and collect end time
% files.QVanalysis = [files.path,'script_QVanalysis__xy-sym.lsf']; % Lumerical script to calculate partial Qs and V
files.QVanalysis = [files.path,'script_QVanalysis__y-sym.lsf']; % Lumerical script to calculate partial Qs and V
% files.QVanalysis = [files.path,'script_QVanalysis__no-sym.lsf']; % Lumerical script to calculate partial Qs and V
files.plotfields = [files.path,'script_plotfields.lsf']; % Lumerical script to plot field profiles

% directory path for data location
datLoc = [pwd,'/trials/EKFreeGeom062519_01/']; 

disp(datLoc)

if ~exist(datLoc,'dir')
    mkdir(datLoc)
end


%% run FDTD simulation
addpath([pwd,'/Nanocavity optimization scripts/']);
ods = runNanobeamCavity(P,files,datLoc);
 
%  heights = [130,150,160,180,200]*10^-9;
%  widths = [130,150,160,180,200]*10^-9;
%  for hhL = heights
%      P.hhL = hhL;
%      for hhR = heights
%          P.hhR = hhR;
%          for hwL = widths
%              P.hwL = hwL;
%              for hwR = widths
%                  P.hwR = hwR;
%                  ods = runNanobeamCavity(P,files,datLoc);
%              end
%          end
%      end
%  end
