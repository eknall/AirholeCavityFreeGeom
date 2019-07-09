% This function runs an FDTD simulation for the TE-like modes of a single nanobeam cavity.
% The input structure P is assumed to have the following entries (plus additional
% geometric parameters, see createNanobeamCavity.m):
%
% M.J. Burek, 08/16

function ds = runNanobeamCavity(P,file,datLoc)
fname = ['PhC-tri_',num2str(P.theta),'o_',...
    'aL_',num2str(P.aL*1e9,'%.0f'),'nm_',...
    'aR_',num2str(P.aR*1e9,'%.0f'),'nm_',...
    'w_',num2str(P.w*1e9,'%.0f'),'nm_', ...
    'hhL_',num2str(P.hhL*1e9,'%.0f'),'nm_', ...
    'hwL_',num2str(P.hwL*1e9,'%.0f'),'nm_', ...
    'hhR_',num2str(P.hhR*1e9,'%.0f'),'nm_', ...
    'hwR_',num2str(P.hwR*1e9,'%.0f'),'nm_', ...
    'nholes_',num2str(P.nholes,'%.0f'),'_', ...
    'ndef_',num2str(P.ndef,'%.0f'),'_', ...
    'maxdef_',num2str(P.maxdef,'%.4f'),'_' ...
    'oblong_',num2str(P.oblong,'%.4f')];

%% Set up geometry
P = createNanobeamCavity(P);
plotdefectcells(P);
pathFig = [datLoc,fname,'_geom.png'];
saveas(gcf,pathFig);
close;

ds.P = P; 

%% Initialize Lumerical FDTD simulation file
FDTDexeLoc = '/n/home08/eknall/opt/lumerical/fdtd/bin/fdtd-solutions';
FDTDmpiLoc = '/n/home08/eknall/opt/lumerical/fdtd/bin/fdtd-engine-impi-lcl';



template = file.template;
[~,templatename,~] = fileparts(template);
build = file.build;
runsim = file.runsim;
QVanalysis = file.QVanalysis;
plotfields = file.plotfields;

%Copy template file to data folder and rename
copyfile(template,datLoc);
movefile([datLoc,templatename,'.fsp'],[datLoc,fname,'.fsp'])

%% Build nanobeam cavity geometry
nbeam = P.nbeam;
w = P.w;
theta = P.theta;
asym = P.asym;
len = P.beamLen;
targLambda = P.lambda;

hh = P.geom(:,1);
hw = P.geom(:,2);
xpos = P.geom(:,3);
ypos = P.geom(:,4);

ahole = P.ahole;

beamStruct=fopen([datLoc,'beamStruct.txt'],'w+');
fprintf(beamStruct,'%1.4f %1.12f %2.1f %1.12f %1.12f %1.12f\r\n', [nbeam;w;theta;len;asym;targLambda]);
fclose(beamStruct);

holeStruct=fopen([datLoc,'holeStruct.txt'],'w+');
numHoles = length(hh);
for k = 1:numHoles
    fprintf(holeStruct,'%1.12f %1.12f %1.12f %1.12f \r\n', [hh(k);hw(k);xpos(k);ypos(k)]);
end
fclose(holeStruct);

periodStruct=fopen([datLoc,'periodStruct.txt'],'w+');
numHoles = length(ahole);
for k = 1:numHoles
    fprintf(periodStruct,'%1.12f \r\n', ahole(k));
end
fclose(periodStruct);

eval(['!"',FDTDexeLoc,'" "',[datLoc,fname,'.fsp'],'" -nw -run "',build,'" -exit'])

% delete([datLoc,'beamStruct.txt']); %erase tmp .txt files
% delete([datLoc,'holeStruct.txt']);

% copy assembled FDTD file and append '_solved.fsp' before running
% simulation
copyfile([datLoc,fname,'.fsp'],[datLoc,fname,'_solved.fsp']);

%% Run FDTD simulation and save
tic;
%run on cluster
eval(['! mpirun -np 50 "',FDTDmpiLoc,'" "',[datLoc,fname,'_solved.fsp'],'" '])
eval(['!"',FDTDexeLoc,'" "',[datLoc,fname,'_solved.fsp'],'" -nw -run "',runsim,'" -exit'])
% eval(['!"',FDTDexeLoc,'" "',[datLoc,fname,'_solved.fsp'],'" -run "',runsim,'" '])
load([datLoc,'realSol.mat']);
delete([datLoc,'realSol.mat']); %erase tmp .mat file
ds.simtime = toc; % seconds (simulation time)

tmp = dir([datLoc,fname,'_solved.fsp']); % collect file size
ds.filesize = tmp.bytes; %megabytes
if realSol == 0
    disp('No nanobeam cavity mode found for this geometry')
    
  %  delete([datLoc,fname,'.fsp']); % delete .fsp file
   % delete([datLoc,fname,'_solved.fsp']); % delete solved .fsp file
  %  delete([datLoc,fname,'_solved_p0.log']); % delete solved .fsp file log
  %  delete([datLoc,fname,'_geom.png']); % delete geometry figure
    
    disp('...deleted bad iteration files')
else
    disp('Nanobeam cavity mode exists')
end
disp(['FDTD simulation took ', num2str(ds.simtime/60,'%.2f'),'min and used ',num2str(ds.filesize*1e-6,'%.3f'), 'MB of diskspace'])

%% Post processing
if realSol == 1
    ds.realSol = 1;
    
    % calculate nanobeam cavity partial Qs & V
    eval(['!"',FDTDexeLoc,'" "',[datLoc,fname,'_solved.fsp'],'" -nw -run "',QVanalysis,'" -exit'])
    QVdat = load([datLoc,'QVanalysis.mat']);
    delete([datLoc,'QVanalysis.mat']); %erase .txt files
    
    % plot nanobeam cavity field profiles
    eval(['!"',FDTDexeLoc,'" "',[datLoc,fname,'_solved.fsp'],'" -nw -run "',plotfields,'" -exit'])
    fielddat = load([datLoc,'fieldprofiles.mat']);
    delete([datLoc,'fieldprofiles.mat']); %erase .txt files
    indx3d = load([datLoc,'index.mat']);
    
    plotfieldprofilesXY(P,fielddat,indx3d,QVdat);
    pathFig = [datLoc,fname,'_modeXY.png'];
    saveas(gcf,pathFig);
    close;
    
    plotfieldprofilesYZ(P,fielddat,QVdat);
    pathFig = [datLoc,fname,'_modeYZ.png'];
    saveas(gcf,pathFig);
    close;
    
    % assemble results
    ds.lambda = QVdat.lambda;
    ds.Vmode = real(QVdat.Vmode);
    ds.Qtime = QVdat.Qtime;
    ds.Qt = QVdat.Qt;
    ds.Qwvg = QVdat.Qwvg;
    ds.Qsc = QVdat.Qsc;
    ds.Trans = QVdat.Trans;
    ds.fielddat = fielddat;
    ds.Qz1 = QVdat.Qz1;
    ds.Qz2 = QVdat.Qz2;
    ds.Qx1 = QVdat.Qx1;
    ds.Qx2 = QVdat.Qx2;
    ds.Qy = QVdat.Qy;
    
    % Save data (remove solved Lumerical .fsp to limit harddisk use)
    if P.storeFDTD == 1
        delete([datLoc,fname,'.fsp']); % delete unsolved .fsp file
    else
        delete([datLoc,fname,'_solved.fsp']); % delete solved .fsp file
        delete([datLoc,fname,'_solved_p0.log']); % delete solved .fsp file log
    end
    % save matlab data
    pathMat = [datLoc,fname,'.mat'];
    save(pathMat,'ds');
    
else
    ds.realSol = [];
end
