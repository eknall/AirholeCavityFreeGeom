%updated 07/09/19

function pp = optimizeNanobeamCavity
global P;
global files;
global datLoc;
%% Global parameters

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
P.ndef = 5;                               % # of holes in 1/2 defect region

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

P.storeFDTD = 0;                        % 1 to keep solved FDTD file


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


fnameBase = ['PhC-Asym_',num2str(P.theta),'o_', ...
    'aL0_',num2str(P.aL*1e9,'%.0f'),'nm_', ...
    'aR0_',num2str(P.aR*1e9,'%.0f'),'nm_', ...
    'nholes_',num2str(P.nholes,'%.0f'),'_', ...
    'ndef_',num2str(P.ndef,'%.0f')];

if P.consthole == 1
    fnameBase = [fnameBase,'_consthole'];
end

% directory path for data location
datLoc = ['D:\Files\FDTD simulations\Optimization runs\scattering design trial\',fnameBase,'\']; %%%%%%%%%%%%%%%%%%%%%%%%%%% Update for location of simulation files
datLoc = [pwd,'/trials/EK062619_Start_280-250-150-150-160-160_00/'];

if ~exist(datLoc,'dir')
    mkdir(datLoc)
end


%% Execute search for optimized design

% Initialize pseudorandom number generator state
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
global itrPath;
global y;

%xmax = [0.65,0.55,0.3*(1-P.consthole),0.2,0.550];
%xmin = [0.35,0.25,0.1*(1-P.consthole),0.1,0.350];
%xnom = [P.hh/P.a,P.hw/P.w,0.1*P.oblong*(1-P.consthole),P.maxdef,P.w*1e6];
xmax = [0.330,0.330,0.78,0.90,0.78,0.90,0.550];
xmin = [0.210,0.210,0.1,0.1,0.1,0.1,0.480];
xnom = [P.aL*1e6,P.aR*1e6,P.hhL/P.aL,P.hwL/P.w,P.hhR/P.aR,P.hwR/P.w,P.w*1e6];
randomStart = 0; % 1 if we want to pick a random starting point Random start not modded to work
% with asymmetric cavity. will crash
strPoint = 1; % starting point tally
finished = 0;

while ~finished
    license = 0;
    disp('test text')
    while ~license
        y = zeros(size(xmin));
        
        if randomStart == 1
            % Pick a random point within the parameter space, but ensure that it
            % has modes around our target before proceeding with the optimization
            goodStart = 0;
            disp('Finding starting point');
            while ~goodStart
                xinit = (0.8*rand(size(xmax))+0.2).*(xmax-xmin)+xmin; %sample on 20-80% of the bounds
                disp(xinit);
                P.hh = xinit(1)*P.a;
                P.hw = xinit(2)*xinit(5)*1e-6;
                P.oblong = xinit(3)*10;
                P.maxdef = xinit(4);
                P.w = xinit(5)*1e-6;
                
                P.th = P.w/(2*tan(P.theta*pi/180));   % beam thickness
                
                fname = ['PhC-tri_',num2str(P.theta),'o_',...
                    'a_',num2str(P.aL*1e9,'%.0f'),'nm_',...
                    'w_',num2str(P.w*1e9,'%.0f'),'nm_', ...
                    'hh_',num2str(P.hh*1e9,'%.0f'),'nm_', ...
                    'hw_',num2str(P.hw*1e9,'%.0f'),'nm_', ...
                    'nholes_',num2str(P.nholes,'%.0f'),'_', ...
                    'ndef_',num2str(P.ndef,'%.0f'),'_', ...
                    'maxdef_',num2str(P.maxdef,'%.4f'),'_', ...
                    'oblong_',num2str(P.oblong,'%.4f')];
                
                try
                    ods = runNanobeamCavity(P,files,datLoc);
                    
                    % Only proceed if we have an optical mode in the
                    % ballpark of our target
                    if ~isempty(ods.realSol)
                        % check initial fitness function, try again if < 0.005
                        % this should avoid 3rd order modes
                        if ((min(ods.Qsc,1e6)/1e6)/ods.Vmode) > 0.005
                            goodStart = 1;
                        else
                            disp('Inadequate starting point, retrying...')
                            delete([datLoc,fname,'.mat']); % delete solved .mat file
                            delete([datLoc,fname,'_geom.png']); % delete geometry figure
                            delete([datLoc,fname,'_modeXY.png']); % delete mode profile figure
                            delete([datLoc,fname,'_modeYZ.png']); % delete mode profile figure
                            if P.storeFDTD == 1 % remove leftover lumerical files
                                delete([datLoc,fname,'_solved.fsp']);
                                delete([datLoc,fname,'_solved_p0.log']);
                            else
                                delete([datLoc,fname,'.fsp']); % delete unsolved .fsp file
                            end
                        end
                    end
                catch lasterror
                    disp(lasterror);
                end
            end
            disp('Starting point found');
            strPoint = strPoint + 1;
        else
            xinit = xnom;
            strPoint = 1;
        end
        
        % create .txt file to assemble iteration results
        itrPath = [datLoc,fnameBase,'_starting point ',num2str(strPoint),'.txt'];
        itr = fopen(itrPath,'wt+');
        fclose(itr);
        
        try
            pp = fminsearchbnd(@(x) -1*fitness(x),xinit,xmin,xmax);
            finished = 1; %fminsearch completes before min step (rare)
            license = 1;
        catch lasterror %will occur when fminsearch exits (min step or license NA)
            disp(lasterror);
            finished = 0;
            if strcmpi(lasterror.message,'minimum step size hit')
                % when min step hits, exit license search
                license = 1;
                %Write file to indicate reached minimum step size
                indicatorFilePath = [datLoc,fnameBase,'_starting point ',num2str(strPoint),'_minStepReached.txt'];
                indicF = fopen(indicatorFilePath,'wt+');
                fclose(indicF);
                % when min step hits, stop optimization if not using random
                % starting points
                if randomStart == 0
                    finished = 1; %avoid multiple optimization passes when not necessary
                end
            else
                disp('License not available')
            end
        end
    end
end
end



function F = fitness(x)
% x(1) = hh (in as)
% x(2) = hw (in widths)
% x(3) = oblong
% x(4) = maxdef
% x(5) = width (in um)
%EK mod for asym
% x(1) = aL (in um)
% x(2) = aR
% x(3) = hhL (in as)
% x(4) = hwL (in widths)
% x(5) = hhR
% x(6) = hwR
% x(7) = width (in um)
global P;
global files;
global datLoc;
global itrPath;
global y;

%P.hh = x(1)*P.a;
%P.hw = x(2)*x(5)*1e-6;
%P.oblong = x(3)*10;
%P.maxdef = x(4);
%P.w = x(5)*1e-6;
%P.th = P.w/(2*tan(P.theta*pi/180));   % beam thickness

P.aL = x(1)*1e-6;
P.aR = x(2)*1e-6;
P.hhL = x(3)*P.aL;
P.hwL = x(4)*x(5)*1e-6;
P.hhR = x(5)*P.aR;
P.hwR = x(6)*x(5)*1e-6;
P.w = x(7)*1e-6;
P.th = P.w/(2*tan(P.theta*pi/180));   % beam thickness



disp(x);
if(max(abs(x-y)./y)<0.001)
    error('minimum step size hit'); %exit fminsearch
end
disp(['geometry change: ',num2str(max(abs(x-y)./y)*100),'%'])
y = x;

try
    disp('Solving for localized optical mode')
    ods = runNanobeamCavity(P,files,datLoc);
    
    % Only proceed if we have an optical mode in the ballpark of our target
    if isempty(ods.realSol)
        disp(['No optical modes near target wavelength ',num2str(P.lambda*1e9,'%.0f'),' nm']);
        F = 0; %continue fminsearch
    else
        % Once we have found our nanobeam cavity mode for this
        % structure, calculate the fitness function
        
        % This ensures we only keep very high Q modes, allowing
        % wvg losses to be engineered once an optimal design is
        % converged upon. Anything above about 500,000 is irrelevant
        % since we won't actually be able to achieve it due to fab
        % errors, surface effects, the phase of the moon, etc.
        Qfac = (min(ods.Qt,1e6));
        Qfacx1 = (min(ods.Qx1,1e6));
        Qfacx2 = (min(ods.Qx2,1e6));
        QfacSc = (min(ods.Qsc,1e6));
        
        
        % In addition to good cooperativity, we want preferential light
        % emission to the waveguide, so we include this in the fitness
        % parameter 
        Qwvgs = 1/(1/Qfacx1 + 1/Qfacx2);
		guidedNess = QfacSc/Qwvgs;
		rightsidedness = Qfacx2/Qfacx1;
		wavelenPenalty = exp(-(5^-2)*(ods.lambda*10^9-P.lambda*10^9)^2);
		F = guidedNess*rightsidedness/ods.Vmode;
		F = F*wavelenPenalty*Qfac*QfacSc/ods.Vmode;
        F = sqrt(F);
        disp(['fitness value = ',num2str(F)]);
        
        % Append iteration fitness parameter in .txt file
        itr = fopen(itrPath,'at+');
        fprintf(itr,['%.0f %.3e %.3e %.0f %.0f', ...
            ' %.6e %.6e %.6e %.6e %.6e', ...
            ' %.4f %.4f', ...
            ' %.6e', ...
            ' %.6e', ...
            ' %.6e', ...
            ' %.6e', ...
            ' %.6e', ...
            ' %.6e', ...
            ' %.6f', ...
            ' %.6e', ...
            ' %.6f', ...
            ' %.6e', ...
            ' %.6f', ...
            ' %.6f\r\n'],...
            P.theta,P.aL,P.aR,P.nholes,P.ndef,...
            P.w,P.hhL,P.hwL,P.hhR,P.hwR,...
            P.maxdef,P.oblong,...
            ods.lambda, ...
            ods.Qtime, ...
            ods.Qt, ... 
            ods.Qsc, ...
            ods.Qz1, ...
            ods.Qz2, ...
            ods.Qx1, ...
            ods.Qx2, ...
            ods.Qy, ...
            ods.Trans, ...
            ods.Vmode, ...
            F);
        fclose(itr);
    end
catch lasterror
    if strcmpi(lasterror.identifier,'MATLAB:load:couldNotReadFile')
        error(lasterror) %exit fminsearch
    end
    disp(lasterror); %display litho constraint
    F = 0; %continue fminsearch
end

end



function [x,fval,exitflag,output]=fminsearchbnd(fun,x0,LB,UB,options,varargin)
% FMINSEARCHBND: FMINSEARCH, but with bound constraints by transformation
% usage: x=FMINSEARCHBND(fun,x0)
% usage: x=FMINSEARCHBND(fun,x0,LB)
% usage: x=FMINSEARCHBND(fun,x0,LB,UB)
% usage: x=FMINSEARCHBND(fun,x0,LB,UB,options)
% usage: x=FMINSEARCHBND(fun,x0,LB,UB,options,p1,p2,...)
% usage: [x,fval,exitflag,output]=FMINSEARCHBND(fun,x0,...)
%
% arguments:
%  fun, x0, options - see the help for FMINSEARCH
%
%  LB - lower bound vector or array, must be the same size as x0
%
%       If no lower bounds exist for one of the variables, then
%       supply -inf for that variable.
%
%       If no lower bounds at all, then LB may be left empty.
%
%       Variables may be fixed in value by setting the corresponding
%       lower and upper bounds to exactly the same value.
%
%  UB - upper bound vector or array, must be the same size as x0
%
%       If no upper bounds exist for one of the variables, then
%       supply +inf for that variable.
%
%       If no upper bounds at all, then UB may be left empty.
%
%       Variables may be fixed in value by setting the corresponding
%       lower and upper bounds to exactly the same value.
%
% Notes:
%
%  If options is supplied, then TolX will apply to the transformed
%  variables. All other FMINSEARCH parameters should be unaffected.
%
%  Variables which are constrained by both a lower and an upper
%  bound will use a sin transformation. Those constrained by
%  only a lower or an upper bound will use a quadratic
%  transformation, and unconstrained variables will be left alone.
%
%  Variables may be fixed by setting their respective bounds equal.
%  In this case, the problem will be reduced in size for FMINSEARCH.
%
%  The bounds are inclusive inequalities, which admit the
%  boundary values themselves, but will not permit ANY function
%  evaluations outside the bounds. These constraints are strictly
%  followed.
%
%  If your problem has an EXCLUSIVE (strict) constraint which will
%  not admit evaluation at the bound itself, then you must provide
%  a slightly offset bound. An example of this is a function which
%  contains the log of one of its parameters. If you constrain the
%  variable to have a lower bound of zero, then FMINSEARCHBND may
%  try to evaluate the function exactly at zero.
%
%
% Example usage:
% rosen = @(x) (1-x(1)).^2 + 105*(x(2)-x(1).^2).^2;
%
% fminsearch(rosen,[3 3])     % unconstrained
% ans =
%    1.0000    1.0000
%
% fminsearchbnd(rosen,[3 3],[2 2],[])     % constrained
% ans =
%    2.0000    4.0000
%
% See test_main.m for other examples of use.
%
%
% See also: fminsearch, fminspleas
%
%
% Author: John D'Errico
% E-mail: woodchips@rochester.rr.com
% Release: 4
% Release date: 7/23/06

% size checks
xsize = size(x0);
x0 = x0(:);
n=length(x0);

if (nargin<3) || isempty(LB)
    LB = repmat(-inf,n,1);
else
    LB = LB(:);
end
if (nargin<4) || isempty(UB)
    UB = repmat(inf,n,1);
else
    UB = UB(:);
end

if (n~=length(LB)) || (n~=length(UB))
    error 'x0 is incompatible in size with either LB or UB.'
end

% set default options if necessary
if (nargin<5) || isempty(options)
    options = optimset('fminsearch');
end

% stuff into a struct to pass around
params.args = varargin;
params.LB = LB;
params.UB = UB;
params.fun = fun;
params.n = n;
params.OutputFcn = [];

% 0 --> unconstrained variable
% 1 --> lower bound only
% 2 --> upper bound only
% 3 --> dual finite bounds
% 4 --> fixed variable
params.BoundClass = zeros(n,1);
for i=1:n
    k = isfinite(LB(i)) + 2*isfinite(UB(i));
    params.BoundClass(i) = k;
    if (k==3) && (LB(i)==UB(i))
        params.BoundClass(i) = 4;
    end
end

% transform starting values into their unconstrained
% surrogates. Check for infeasible starting guesses.
x0u = x0;
k=1;
for i = 1:n
    switch params.BoundClass(i)
        case 1
            % lower bound only
            if x0(i)<=LB(i)
                % infeasible starting value. Use bound.
                x0u(k) = 0;
            else
                x0u(k) = sqrt(x0(i) - LB(i));
            end
            
            % increment k
            k=k+1;
        case 2
            % upper bound only
            if x0(i)>=UB(i)
                % infeasible starting value. use bound.
                x0u(k) = 0;
            else
                x0u(k) = sqrt(UB(i) - x0(i));
            end
            
            % increment k
            k=k+1;
        case 3
            % lower and upper bounds
            if x0(i)<=LB(i)
                % infeasible starting value
                x0u(k) = -pi/2;
            elseif x0(i)>=UB(i)
                % infeasible starting value
                x0u(k) = pi/2;
            else
                x0u(k) = 2*(x0(i) - LB(i))/(UB(i)-LB(i)) - 1;
                % shift by 2*pi to avoid problems at zero in fminsearch
                % otherwise, the initial simplex is vanishingly small
                x0u(k) = 2*pi+asin(max(-1,min(1,x0u(k))));
            end
            
            % increment k
            k=k+1;
        case 0
            % unconstrained variable. x0u(i) is set.
            x0u(k) = x0(i);
            
            % increment k
            k=k+1;
        case 4
            % fixed variable. drop it before fminsearch sees it.
            % k is not incremented for this variable.
    end
    
end
% if any of the unknowns were fixed, then we need to shorten
% x0u now.
if k<=n
    x0u(k:n) = [];
end

% were all the variables fixed?
if isempty(x0u)
    % All variables were fixed. quit immediately, setting the
    % appropriate parameters, then return.
    
    % undo the variable transformations into the original space
    x = xtransform(x0u,params);
    
    % final reshape
    x = reshape(x,xsize);
    
    % stuff fval with the final value
    fval = feval(params.fun,x,params.args{:});
    
    % fminsearchbnd was not called
    exitflag = 0;
    
    output.iterations = 0;
    output.funcount = 1;
    output.algorithm = 'fminsearch';
    output.message = 'All variables were held fixed by the applied bounds';
    
    % return with no call at all to fminsearch
    return
end

% Check for an outputfcn. If there is any, then substitute my
% own wrapper function.
if ~isempty(options.OutputFcn)
    params.OutputFcn = options.OutputFcn;
    options.OutputFcn = @outfun_wrapper;
end

% now we can call fminsearch, but with our own
% intra-objective function.
[xu,fval,exitflag,output] = fminsearch(@intrafun,x0u,options,params);

% undo the variable transformations into the original space
x = xtransform(xu,params);

% final reshape
x = reshape(x,xsize);

% Use a nested function as the OutputFcn wrapper
    function stop = outfun_wrapper(x,varargin);
        % we need to transform x first
        xtrans = xtransform(x,params);
        
        % then call the user supplied OutputFcn
        stop = params.OutputFcn(xtrans,varargin{1:(end-1)});
        
    end

end % mainline end

% ======================================
% ========= begin subfunctions =========
% ======================================
function fval = intrafun(x,params)
% transform variables, then call original function

% transform
xtrans = xtransform(x,params);

% and call fun
fval = feval(params.fun,xtrans,params.args{:});

end % sub function intrafun end

% ======================================
function xtrans = xtransform(x,params)
% converts unconstrained variables into their original domains

xtrans = zeros(1,params.n);
% k allows some variables to be fixed, thus dropped from the
% optimization.
k=1;
for i = 1:params.n
    switch params.BoundClass(i)
        case 1
            % lower bound only
            xtrans(i) = params.LB(i) + x(k).^2;
            
            k=k+1;
        case 2
            % upper bound only
            xtrans(i) = params.UB(i) - x(k).^2;
            
            k=k+1;
        case 3
            % lower and upper bounds
            xtrans(i) = (sin(x(k))+1)/2;
            xtrans(i) = xtrans(i)*(params.UB(i) - params.LB(i)) + params.LB(i);
            % just in case of any floating point problems
            xtrans(i) = max(params.LB(i),min(params.UB(i),xtrans(i)));
            
            k=k+1;
        case 4
            % fixed variable, bounds are equal, set it at either bound
            xtrans(i) = params.LB(i);
        case 0
            % unconstrained variable.
            xtrans(i) = x(k);
            
            k=k+1;
    end
end

end % sub function xtransform end





