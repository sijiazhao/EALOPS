clear;close all;

% experimentname = 'Expt1 - Young'; % Uncomment to plot figure 3
experimentname = 'Expt2 - Old'; % Uncomment to plot figure 7

pupilunit = 'mm';

path_root = 'E:\Sijia\EALOPs\Submission\Upload\Processed data_Main';
path_in = [path_root '\' experimentname];

blocklist = 1:6;
trial_number = 12;
condlist = {'Easy','Medium','Hard'};%{'1 stream','2 streams','3 streams'};
%for plotting PDR we realign
tw_epoch = [-2 29]; % time window for epoch [s]
tw_bc = [3.5,4];
tw_PDR = [-6,25];
tw_bc_PDR = [-0.5,0];
tw_stats = [-0.5,25];

%% Retrieve subject informtion from excel file
[~,~,expMat] = xlsread([path_in '\subject_Info.xlsx']);
FsMatrix = expMat(:,4); %sampling frequency, you can set it as 1000Hz
sublist = expMat(:,1); %cell array for subject names
smpfreq = FsMatrix{1};
timeaxis = (tw_epoch(1)*smpfreq:tw_epoch(2)*smpfreq)/smpfreq;
define_colourmap;

filename = [path_in, '\behav_result.mat'];
load(filename);

filename = [path_in, '\sublist_group.mat'];
load(filename,'sublist_good','sublist_poor','sublist_exclude','trial2keep_valid');
sublist_include = [sublist_good;sublist_poor];

%% Prepare the time-binned behavioural result
bin_label = {'0-5s','5-10s','10-15s','15-20s','20-25s'};
timeaxis_bin = 1:numel(bin_label);

edges = (0:5:25);
edges1 = edges+tw_bc(2);
data = [];
for s = 1:numel(sublist)
    for k = 1:numel(condlist)
        bot = over_time_Hits_Target{s,k};
        for b = 1:(numel(edges1)-1)
            [~,tw(1)] = fFindClosest(timeaxis,edges1(b));
            [~,tw(2)] = fFindClosest(timeaxis,edges1(b+1));
            tw = tw(1):tw(2)-1;
            m = bot(:,tw);
            data(s,k,b) = (numel(find(m==1)))/(numel(find(m==1))+ numel(find(m==-1)))*100;
        end
    end
end
BOT_hit = data;
BOT_hit(sublist_exclude,:,:) = NaN;

data = [];
for s = 1:numel(sublist)
    for k = 1:numel(condlist)
        bot = over_time_Fas{s,k};
        for b = 1:(numel(edges1)-1)
            [~,tw(1)] = fFindClosest(timeaxis,edges1(b));
            [~,tw(2)] = fFindClosest(timeaxis,edges1(b+1));
            tw = tw(1):tw(2)-1;
            m = bot(:,tw);
            data(s,k,b) = numel(find(m == 1));
        end
    end
end
BOT_fa = data;
BOT_fa(sublist_exclude,:,:) = NaN;

%% Start to plot
figure(1); clf;

% Time-binned result: basic hit rate
subplot(7,3,[7,8]); hold on;
hold on;
for cond = 1:numel(condlist)
    cond_mean = nanmean(squeeze(BOT_hit(:,cond,:)),1);
    cond_sem = nanstd(squeeze(BOT_hit(:,cond,:)),1)/sqrt(length(sublist_include));
    
    errorbar(timeaxis_bin,cond_mean,cond_sem,'Marker','none','Color',colourmap(cond,:),'LineWidth',2,'CapSize',10);
end

xlim([0.75 5.25]);
ylim([70 100]);
ax=gca;
ax.XTick = 1:1:numel(bin_label);
ax.XTickLabel = bin_label;
set(gca,'XTick',[]);
set(gca, 'YTick', 70:5:100);
ylabel('Hit rate [%]');
hold off;

% Time-binned result: false alarm rate
subplot(7,3,[10,11]); hold on;
for cond = 1:numel(condlist)
    cond_mean = nanmean(squeeze(BOT_fa(:,cond,:)),1);
    cond_sem = nanstd(squeeze(BOT_fa(:,cond,:)),1)/sqrt(length(sublist_include));
    
    errorbar(timeaxis_bin,cond_mean,cond_sem,'Marker','none','Color',colourmap(cond,:),'LineWidth',2,'CapSize',10);
end
xlim([0.75 5.25]);
ylim([0 5]);
ax=gca;
ax.XTick = 1:1:numel(bin_label);
ax.XTickLabel = bin_label;
ax.YTick = 0:1:5;
ylabel('#False alarm');
hold off;

% Hit rate difference between hard and medium (time-binned)
diff_hit_2vs3 = squeeze(BOT_hit(:,3,:)) - squeeze(BOT_hit(:,2,:));
subplot(7,3,[13,14]); hold on;

for t = 1:numel(timeaxis_bin)
    y=diff_hit_2vs3(:,t);
    scatter(timeaxis_bin(t)*ones(size(y)),y,'filled','MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor','none','LineWidth',0.5);
end
cond_mean = nanmean(diff_hit_2vs3,1);
% cond_sem = nanstd(diff_hit_2vs3,1)/sqrt(length(sublist_include));
cond_sem = nanstd(diff_hit_2vs3,1);
errorbar(timeaxis_bin,cond_mean,cond_sem,'-o','Color',colourmap_group(2,:),'LineWidth',2,'CapSize',10);

xlim([0.75 5.25]);
ylim([-60 40]);
ax=gca;
ax.XTick = 1:1:numel(bin_label);
ax.XTickLabel = bin_label;
set(gca,'XTick',[]);
set(gca, 'YTick', -60:20:40);
ylabel({'Diff. in Hit rate'; 'Hard-Medium [%]'});
hold off;


%% Prepare PDR result
file_in = [path_in,'\pupil_PDRatOnset_(nor0,bc1)hir.mat'];
load(file_in);

timeaxis = (tw_PDR(1)*1000:tw_PDR(2)*1000)/1000;

% Z-score
tw_bc0 = [3.5,4];
idx_bc0 = find(timeaxis>=(tw_bc0(1)) & timeaxis<=(tw_bc0(2)));

tw_bc1 = [-5,0];
idx_bc1 = find(timeaxis>=(tw_bc1(1)) & timeaxis<=(tw_bc1(2)));

if strcmp(pupilunit(1),'z')
    for subj = 1:numel(sublist)
        pp = [];
        for cond = 1:numel(condlist)
            p = P{subj,cond};
            pp = [pp;p];
        end
        
        for cond = 1:numel(condlist)
            p = P{subj,cond};
            switch pupilunit
                case 'z(full)'
                    bas = pp;
                case 'z(baseline)'
                    bas = pp(:,idx_bc0);
                case 'z(longbaseline)'
                    bas = pp(:,idx_bc1);
            end
            
            bas = reshape(bas,[numel(bas),1]);
            
            
            p = (p-nanmean(bas))/nanstd(bas); % Normalised to z score
            P{subj,cond} = p;
        end
    end
    disp('Pupil data normalised');
end

% Compute PDR
pupilmean = [];
for subj = 1:length(sublist)
    for cond = 1:numel(condlist)
        x = P{subj,cond};
        pupilmean(subj,cond,:) = nanmean(x,1);
    end
end

% Convert unit
switch pupilunit
    case 'mm'
        % diameter of the black dot = 5mm
        dot.left = 795.6;
        dot.right = 766.2;
        dot.both = 793;
        pupilmean = pupilmean/dot.left;
end

% Smooth again for plotting
for subj = 1:length(sublist)
    for cond = 1:numel(condlist)
        Pn = pupilmean(subj,cond,:);
        Pn = smooth(Pn,150,'hann'); %Smooth pupil data
        pupilmean(subj,cond,:) = Pn;
    end
end

% Downsample
fsnew = 20;
pupilmean = run4ds(pupilmean,timeaxis,fsnew);
timeaxis = tw_PDR(1):1/fsnew:tw_PDR(2);

subplot(7,3,[1,2,4,5]); hold on;
for k = 1:numel(condlist)
    indmean = squeeze(pupilmean(:,k,:));
    a = nanmean(indmean);
    plot(timeaxis,a,'LineWidth',1.5,'Color',colourmap(k,:));
    
    % PLOT error bar
    a = nanmean(indmean)';
    b = (nanstd(indmean)/sqrt(size(indmean,1)))';
    
    curve1 = a+b;
    curve2 = flipud(a-b);
    
    X = [timeaxis'; flipud(timeaxis')];
    Y = [curve1; curve2];
    figfill = fill(X,Y,colourmap(k,:),'edgecolor','none','facealpha',0.2);
    
    set(gca,'children',flipud(get(gca,'children'))); % Send shaded areas in background
end

switch pupilunit
    case 'mm'
        ylim([-.5 .5]);
    case 'z(longbaseline)'
        ylim([-1.2 1]);
end

% PLOT BOOTSTRAP STATS
ab = ylim;
aa= 0.025*(max(ab)-min(ab));
dist = 0.09*(max(ab)-min(ab));
Ystat = ab(1);
for kk = 1:3
    statpairs = [1,2;1,3;2,3];
    
    Ystat = Ystat+aa;
    
    statstime_start = 0;
    statstime_end = max(timeaxis);
    statstime = find((timeaxis>statstime_start) & (timeaxis<=statstime_end));
    cond1 = squeeze(pupilmean(:,statpairs(kk,1),min(statstime):max(statstime)));
    cond2 = squeeze(pupilmean(:,statpairs(kk,2),min(statstime):max(statstime)));
    dataB = bootstrap(cond1-cond2); %Repeated measure
    s_boot = findSigDiff(dataB, 0.05);
    pre_s = NaN(size(find((timeaxis<=statstime_start))));
    s_boot = [pre_s,s_boot];
    statlimit = FindStatLimit(tw_bc_PDR,s_boot,timeaxis);
    statlimit = max(statlimit);
    ClusterBoot;
    s_boot = s_boot_corrected;
    plot(timeaxis,Ystat*abs(s_boot),'LineWidth',5,'Color',colourmap_statspair(kk,:));
end
xlim([-.5 25]);
ylabel(['[' pupilunit ']']);

%% Prepare time-binned PDR result
statspair = [1,2; 1,3; 2,3];
p = squeeze(pupilmean(:,statspair(3,2),:)-pupilmean(:,statspair(3,1),:));
p(sublist_exclude,:) = [];

edges = 0:5:25;
np = [];
for s=1:size(p,1)
    for b = 1:(numel(edges)-1)
        [~,tw(1)] = fFindClosest(timeaxis,edges(b));
        [~,tw(2)] = fFindClosest(timeaxis,edges(b+1));
        tw = tw(1):tw(2);
        np(s,b) = nanmean(p(s,tw));
    end
end
pbin = np;

%% Plot time-binned PDR result
subplot(7,3,[16,17]); hold on;
whichmode = 'time-binned';
y = diff_hit_2vs3;
y(sublist_exclude,:) = [];

RHO = [];
PVAL = [];
for t = 1:numel(timeaxis_bin)
    x = pbin(:,t);
    y1 = y(:,t);
    [RHO(t),PVAL(t)] = corr(x,y1,'Type','Spearman');
end

if_holm_bonferroni = 1;
if if_holm_bonferroni
    sigpoints = test_holm_bonferroni(PVAL);
    
    bar(timeaxis_bin,sigpoints,'r','FaceAlpha',0.5,'EdgeColor','none');
    bar(timeaxis_bin,-sigpoints,'r','FaceAlpha',0.5,'EdgeColor','none');
    %         PVAL(PVAL>(0.05/numel(timeaxis))) = NaN;
else
    PVAL(PVAL>(0.05)) = NaN;
    
    sigpoints = PVAL;
    sigpoints(~isnan(PVAL))= 1;
    sigpoints(isnan(PVAL))= 0;
    bar(timeaxis_bin,sigpoints,'r','FaceAlpha',0.5,'EdgeColor','none');
    sigpoints(~isnan(PVAL))= -1;
    sigpoints(isnan(PVAL))= 0;
    bar(timeaxis_bin,sigpoints,'r','FaceAlpha',0.5,'EdgeColor','none');
end

bar(timeaxis_bin,RHO,'black');

xl = xlim;
line(xl, -0.5*ones(size(xl)),'Color',[.5 .5 .5],'LineStyle','--');
line(xl, 0.5*ones(size(xl)),'Color',[.5 .5 .5],'LineStyle','--');

for t = 1:numel(timeaxis_bin)
    text(timeaxis_bin(t)-.25,.75,sprintf('p=%.4f',PVAL(t)),'FontSize',8);
end
ylim([-1,1]);
ylabel('Correlation coefficient');

xlim([0.75 5.25]);
ax=gca;
ax.XTick = 1:1:numel(bin_label);
ax.XTickLabel = bin_label;
set(gca,'XTick',[]);
set(gca, 'YTick', -1:.5:1);
hold off;


%% Corr PDR with overall behaviour
% load behav results
filename = [path_in,'\behav_result.mat'];
load(filename,'behav');
% filename = [path_in,'\distractScore'];
% load(filename,'distractScore');

statspair = [1,2; 1,3; 2,3];

for subj = 1:length(sublist)
    for k =1:size(statspair,1)
        %hit difference
        data = behav.pHit{1,subj}*100;
        behavcorr_HIT(subj,k) = nanmean(data(:,statspair(k,2)))-nanmean(data(:,statspair(k,1)));
        %FA difference
        data = behav.nFa{1,subj};
        behavcorr_FA(subj,k) = sum(data(:,statspair(k,2)))-sum(data(:,statspair(k,1)));
        %distFA difference
        data = behav.distFa{1,subj};
        behavcorr_distFA(subj,k) = sum(data(:,statspair(k,2)))-sum(data(:,statspair(k,1)));
        %FA rate difference
        data = behav.FaRate{1,subj};
        behavcorr_FARate(subj,k) = nanmean(data(:,statspair(k,2)))-nanmean(data(:,statspair(k,1)));
    end
    data = behav.nFa{1,subj};
    Fa_2and3(subj,1) = sum(data(:,2))+sum(data(:,3));
end

%% Plot the correlation
subplot(7,3,[19,20]); hold on;

whichmode = 'overall';
y = behavcorr_HIT(:,3);
y(sublist_exclude,:) = [];

RHO = [];
PVAL = [];
for t = 1:numel(timeaxis)
    x = p(:,t);
    y = y;
    [RHO(t),PVAL(t)] = corr(x,y,'Type','Spearman');
end


if_holm_bonferroni = 0;
if if_holm_bonferroni
    sigpoints = test_holm_bonferroni(PVAL);
    
    bar(timeaxis,sigpoints,'r','FaceAlpha',0.5,'EdgeColor','none');
    bar(timeaxis,-sigpoints,'r','FaceAlpha',0.5,'EdgeColor','none');
    %         PVAL(PVAL>(0.05/numel(timeaxis))) = NaN;
else
    PVAL(PVAL>(0.05)) = NaN;
    
    sigpoints = PVAL;
    sigpoints(~isnan(PVAL))= 1;
    sigpoints(isnan(PVAL))= 0;
    bar(timeaxis,sigpoints,'r','FaceAlpha',0.5,'EdgeColor','none');
    sigpoints(~isnan(PVAL))= -1;
    sigpoints(isnan(PVAL))= 0;
    bar(timeaxis,sigpoints,'r','FaceAlpha',0.5,'EdgeColor','none');
end

bar(timeaxis,RHO,'black');

xlim([-.5 25]);
xl = xlim;
line(xl, -0.5*ones(size(xl)),'Color',[.5 .5 .5],'LineStyle','--')
line(xl, 0.5*ones(size(xl)),'Color',[.5 .5 .5],'LineStyle','--')


ylim([-1,1]);
set(gca, 'YTick', -1:.5:1);
ylabel('Correlation coefficient','FontSize',10);
xlabel('Time from onset [s]','FontSize',10);
hold off;

%%
subplot(7,3,18); hold on;
t = 4; % the 4th bin
y = diff_hit_2vs3;
y(sublist_exclude,:) = [];
y1 = y(:,t);

x = pbin(:,t);

hold on;
x = x;
y = y1;
scatter(x,y,'o','filled','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','LineWidth', .5);
p = polyfit(x,y,1);
f = polyval(p,x);
plot(x,f,'k');
hold off;

[r,p] = corr(x,y,'Type','Spearman');
xl = [-.5,1.5];
yl = [-40,10];
ylim(yl);
xlim(xl);
text(xl(1),9,['r=',num2str(round(r,4)),',p=',num2str(round(p,4))],'FontSize',8);
xlabel('Diff. in PDR [mm]','FontSize',10);
ylabel('Diff. in HR [%]','FontSize',10);

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 12]);

