clear all; close all; clc
%% Assemble optimization results
datLoc = 'D:\Files\FDTD simulations\Optimization runs\reoptimized designs\';
datLoc = [pwd,'/trials/'];
folders = dir(datLoc);
folders = folders([folders.isdir] & ~strncmpi('.', {folders.name}, 1));
data = {};
for i = 1:length(folders)
    tmp = dir([datLoc,folders(i).name,'/*.txt']);
    for j = 1:length(tmp)
        data{end+1} = [datLoc,folders(i).name,'/',tmp(j).name];
    end
end

%% Revert fitness value definition
Vplot = 0; % set to 1 in order to only plot F = 1/V


%% Data filtering

% filters
hhMin = 100e-9;
hwMin = 100e-9;
Qlow = 0e5;
Vhigh = 2;

results = {};
for i = 1:length(data)
    tmp = load(data{1,i});
    
    % filter small holes
    tmp2 = [];
    for j = 1:length(tmp(:,1))
        if tmp(j,6) > hhMin && ...
                tmp(j,7) > hwMin && ...
                tmp(j,end-2) > Qlow && ...
                tmp(j,end-1) < Vhigh
            tmp2 = [tmp2;tmp(j,:)];
        end
    end
    if ~isempty(tmp2)
        if ~isempty(tmp2(:,1))
            results{end+1} = tmp2;
        end
    end
end

totItr = 0;
maxItr = 0;
maxV = 0;
minV = Inf;
maxF = 0;
minF = Inf;
maxInd = zeros(1,2);
for i = 1:length(results)
    tmp = results{i};
    totItr = totItr + length(tmp(:,1));
    if length(tmp(:,1)) > maxItr
        maxItr = length(tmp(:,1));
    end
    if max(tmp(:,end-1)) > maxV
        maxV = max(tmp(:,end-1));
    end
    if min(tmp(:,end-1)) < minV
        minV = min(tmp(:,end-1));
        if Vplot == 1
            maxInd = [i,find(tmp(:,end-1) == minV)];
        end
    end
    if max(tmp(:,end)) > maxF
        maxF = max(tmp(:,end));
        if Vplot == 0
            maxInd = [i,find(tmp(:,end) == maxF)];
        end
    end
    if min(tmp(:,end)) < minF
        minF = min(tmp(:,end));
    end
end


%% Plot F vs. iteration # for each starting point

% normalized plots
figure; hold on;
cmap = hsv(length(results))./1.5;
for i = 1:length(results)
    tmp = results{i};
    plot(1:length(tmp(:,1)),(i-1)+(tmp(:,end)/max(tmp(:,end))),...
        ':o','markersize',8,'linewidth',2,'Color',cmap(i,:))
    plot(1,(i-1)+(tmp(1,end)/max(tmp(:,end))),...
        'o','markersize',12,'markeredgecolor','r','linewidth',2)
    plot(1:maxItr+1,(i-1)+ones(1,maxItr+1),'--k','linewidth',2);
end
xlim([1 maxItr+1])
ylim([0 0.05+length(results)])
box on

xlabel('# of iterations','fontsize',16)
ylabel('F/max(F)','fontsize',16)
title('F = [min(Q_{sc},Q_{cutoff})/Q_{cutoff}]/V_{mode}','fontsize',16)
set(gca,'fontsize',14,'ytick',[0 1],'yticklabel',{'0' '1'})

% unnormalized plot
figure; hold on;
cmap = hsv(length(results))./1.5;
for i = 1:length(results)
    tmp = results{i};
    plot(1:length(tmp(:,1)),tmp(:,end),...
        ':o','markersize',8,'linewidth',2,'Color',cmap(i,:))
    plot(1,tmp(1,end),...
        'o','markersize',12,'markeredgecolor','r','linewidth',2)
end
ylim([0 1.05*maxF])
xlim([1 maxItr+1])
box on

xlabel('# of iterations','fontsize',16)
ylabel('F','fontsize',16)
title('F = [min(Q_{sc},Q_{cutoff})/Q_{cutoff}]/V_{mode}','fontsize',16)
set(gca,'fontsize',14)


%% Plot F vs. dimension for each starting point
figure;
set(gcf,'position',[969 149 944 1068-200])
axes('position',[0.125 0.1 0.8 0.8])
hold on;

% collect all iteration data
xDatall = [];
yDatall = [];
defDatall = [];
Vall = [];
Fall = [];
Tall =[];
for i = 1:length(results)
    tmp = results{i};
    xDatall = [xDatall; tmp(:,6)./(tmp(:,2))];
    yDatall = [yDatall; tmp(:,7)./tmp(:,5)];
    defDatall = [defDatall; tmp(:,8)];
    Fall = [Fall; tmp(:,end)];
    Vall = [Vall; tmp(:,end-1)];
    Tall = [Tall; ((tmp(:,end-2) - tmp(:,end-3)).^2)./tmp(:,end-2).^2];
end
% sort by ascending F values
[Fall,sort_idx] = sort(Fall);
xDatall = xDatall(sort_idx);
yDatall = yDatall(sort_idx);
defDatall = defDatall(sort_idx);
Vall = Vall(sort_idx);
Tall = Tall(sort_idx);

% change fitness values if only interested in 1/Vmode
if Vplot == 1
    Fall = 1./Vall;
    maxF = 1/minV;
    minF = 1/maxV;
end

% plot shadow points
scatter3(xDatall,yDatall,0.07*ones(1,totItr),10,...
    'markeredgecolor',0.85*[1 1 1],'markerfacecolor',0.85*[1 1 1])
scatter3(0.82*ones(1,totItr),yDatall,defDatall,10,...
    'markeredgecolor',0.85*[1 1 1],'markerfacecolor',0.85*[1 1 1])
scatter3(xDatall,0.82*ones(1,totItr),defDatall,10,...
    'markeredgecolor',0.85*[1 1 1],'markerfacecolor',0.85*[1 1 1])

% plot dashed progress lines
for i = 1:length(results)
    tmp = results{i};
    xDat = tmp(:,6)./tmp(:,2); % hx/a
    yDat = tmp(:,7)./tmp(:,5); % hy/w
    def = tmp(:,8); %maxdef
    
    plot3(xDat,yDat,def,':b',...
        'linewidth',0.5);
    
    % indicate random starting point
    scatter3(xDat(1),yDat(1),def(1),160,'s',...
        'markeredgecolor','r','linewidth',2)
    scatter3(0.82,yDat(1),def(1),100,'s',...
        'markeredgecolor',0.85*[1 1 1],'linewidth',1)
    scatter3(xDat(1),0.82,def(1),100,'s',...
        'markeredgecolor',0.85*[1 1 1],'linewidth',1)
    scatter3(xDat(1),yDat(1),0.07,100,'s',...
        'markeredgecolor',0.85*[1 1 1],'linewidth',1)
end

% plot iterations colored by fitness value
cmap = jet(101);
ind = roundn(0:1e-2:1,-2);
for j = 1:length(Fall)
    for i = 1:length(cmap)
        if roundn((Fall(j)-minF)/(maxF-minF),-2) == ind(i)
            colorP = cmap(i,:);
        end
    end
    scatter3(xDatall(j),yDatall(j),defDatall(j),50,...
        'markeredgecolor',colorP,'markerfacecolor',colorP)
    if sort_idx(j) ~= 1
    text(xDatall(j),yDatall(j),defDatall(j),['   ',num2str(sort_idx(j))])
    end
    
    %indicate global maximum F value
    if Fall(j) == maxF
        scatter3(xDatall(j),yDatall(j),defDatall(j),160,'o',...
            'markeredgecolor','k','linewidth',2)
        scatter3(xDatall(j),yDatall(j),0.07,100,'o',...
            'markeredgecolor',0.75*[1 1 1],'linewidth',1)
        scatter3(0.82,yDatall(j),defDatall(j),100,'o',...
            'markeredgecolor',0.75*[1 1 1],'linewidth',1)
        scatter3(xDatall(j),0.82,defDatall(j),100,'o',...
            'markeredgecolor',0.75*[1 1 1],'linewidth',1)
    end
end

% indicate last iteration
scatter3(xDat(end),yDat(end),tmp(end,8),160,'d',...
    'markeredgecolor','m','linewidth',2)
scatter3(0.82,yDat(end),tmp(end,8),100,'d',...
    'markeredgecolor',0.75*[1 1 1],'linewidth',1)
scatter3(xDat(end),0.82,tmp(end,8),100,'d',...
    'markeredgecolor',0.75*[1 1 1],'linewidth',1)
scatter3(xDat(end),yDat(end),0.07,100,'d',...
    'markeredgecolor',0.75*[1 1 1],'linewidth',1)

axis tight
xlim([0.18 0.82])
ylim([0.18 0.82])
zlim([0.07 0.26])
box on
grid on
view(3)

xlabel('h_x/a','fontsize',20)
ylabel('h_y/w','fontsize',20)
zlabel('d_{max}','fontsize',20)
set(gca,'fontsize',14)

str1 = ['min(F) = ',num2str(minF,'%.2f')];
str2 = ['                                                   max(F) = ',...
    num2str(maxF,'%.2f')];
str3 = [num2str(totItr),' total iterations'];
str4 = [num2str(length(results)),' random starting points'];

cbar = colorbar('northoutside');
set(cbar,'xtick',[minF maxF],'xticklabel', ...
    {str1 str2},...
    'fontsize',12)
cbarPos = get(cbar,'position');
set(cbar,'position',[cbarPos(1)+0.0625 cbarPos(2)-0.0125 0.175 0.015])

if Vplot == 1
    t = title({'F = 1/V_{mode}';[str3,', ',str4]},'fontsize',16);
else
    t = title({'F = [min(Q_{sc},Q_{cutoff})/Q_{cutoff}]/V_{mode}';[str3,', ',str4]},'fontsize',16);
end
tPos = get(t,'position');

%% Collect best iteration

P = results{maxInd(1)}(maxInd(2),:);

fname = ['PhC-tri_',num2str(P(1)),'o_',...
    'a_',num2str(P(2)*1e9,'%.0f'),'nm_',...
    'w_',num2str(P(5)*1e9,'%.0f'),'nm_', ...
    'hh_',num2str(P(6)*1e9,'%.0f'),'nm_', ...
    'hw_',num2str(P(7)*1e9,'%.0f'),'nm_', ...
    'nholes_',num2str(P(3),'%.0f'),'_', ...
    'ndef_',num2str(P(4),'%.0f'),'_', ...
    'maxdef_',num2str(P(8),'%.4f'),'_' ...
    'oblong_',num2str(P(9),'%.4f')];

finLoc = [datLoc,'/optimized design'];

for i = 1:length(folders)
    if exist([datLoc,folders(i).name,'/',fname,'.fsp'],'file')
        if ~exist(finLoc,'dir')
            mkdir(finLoc)
        end
        
        try
            copyfile([datLoc,folders(i).name,'/',fname,'.fsp'],[finLoc,'/',fname,'.fsp'],'f')
            copyfile([datLoc,folders(i).name,'/',fname,'.mat'],[finLoc,'/',fname,'.mat'],'f')
            copyfile([datLoc,folders(i).name,'/',fname,'_geom.png'],[finLoc,'/',fname,'_geom.png'],'f')
            copyfile([datLoc,folders(i).name,'/',fname,'_modeXY.png'],[finLoc,'/',fname,'_modeXY.png'],'f')
            copyfile([datLoc,folders(i).name,'/',fname,'_modeYZ.png'],[finLoc,'/',fname,'_modeYZ.png'],'f')
        catch error
            disp('files already exist')
        end
    end
end


