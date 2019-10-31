clear; close all;

% Run this to plot figure 4
path_in = 'E:\Sijia\EALOPs\Submission\Upload\Processed data_Pupilmetrics';

figure(1);

green = [23,190,37]/255;
purple = [157,23,190]/255;
colourmap_group = [green;purple];

%% Resting pupil diameter measurements
pupilunit = 'mm'; divideByResting = 1;
filename_excel = 'Results_pupil_metrics.xlsx';
load([path_in '\loc_Elder' '[' pupilunit ']'],'pupildiameter1','pupildiameter_std1','pupildiameter2','pupildiameter_std2','PDR','timeaxis','sublist');
p_e_pre = pupildiameter1(:,1)';
p_e_post = pupildiameter2(:,1)';
s_e_pre = pupildiameter_std1(:,1)';
s_e_post= pupildiameter_std2(:,1)';
pdr_e = PDR;

load([path_in '\loc_Young' '[' pupilunit ']'],'pupildiameter1','pupildiameter_std1','pupildiameter2','pupildiameter_std2','PDR','timeaxis','sublist');
p_y_pre = pupildiameter1(:,1)';
p_y_post = pupildiameter2(:,1)';
s_y_pre = pupildiameter_std1(:,1)';
s_y_post= pupildiameter_std2(:,1)';
pdr_y = PDR;

%% ---A---
subplot(5,3,1);
data = [p_y_pre',p_y_post',p_e_pre',p_e_post'];

xlRange = 'A4:D21';
xlswrite(filename_excel,data,xlRange);

plotbar(data,1);
ax=gca;
ax.XTick = 1:1:size(data,2);
ax.XTickLabel = {'Young pre','Young post','Elder pre','Elder post'};
xtickangle(45);
% set(gca, 'YTick', 0:20:100);
ylabel('Pupil diameter [mm]');
% title('Pupil diameter at resting');
% set(gca,'FontSize',12);
xl = [2 10]; %[mm]
ylim(xl);
for x3 = 1:size(data,2)
    p = nanmean(data(:,x3));
    y3 = xl(1)+xl(2)*0.1;%-0.1*(yl(2)-yl(1));
    txt3 = [sprintf('%0.2f',p)];
    text(x3,y3,txt3,'FontSize',9,'HorizontalAlignment','center','Rotation',90);
end

a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Arial','fontsize',8);

%% ---B---
subplot(5,3,2);
data = [s_y_pre',s_y_post',s_e_pre',s_e_post'];

xlRange = 'F4:I21';
xlswrite(filename_excel,data,xlRange);

plotbar(data,1);
ax=gca;
ax.XTick = 1:1:size(data,2);
ax.XTickLabel = {'Young pre','Young post','Elder pre','Elder post'};
xtickangle(45);
ylabel('Variability [mm]');
ylim([0,1.5]);
yl = ylim;
for x3 = 1:size(data,2)
    p = nanmean(data(:,x3));
    y3 = yl(1)+yl(2)*0.1;%-0.1*(yl(2)-yl(1));
    txt3 = [sprintf('%0.2f',p)];
    text(x3,y3,txt3,'FontSize',9,'HorizontalAlignment','center','Rotation',90);
end

a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Arial','fontsize',8);

%% ---C---
subplot(5,3,3);
sn_e_pre = s_e_pre./p_e_pre;
sn_e_post= s_e_post./p_e_post;
sn_y_pre = s_y_pre./p_y_pre;
sn_y_post= s_y_post./p_y_post;

data = [sn_y_pre',sn_y_post',sn_e_pre',sn_e_post'];

xlRange = 'K4:N21';
xlswrite(filename_excel,data,xlRange);

plotbar(data,1);
ax=gca;
ax.XTick = 1:1:size(data,2);
ax.XTickLabel = {'Young pre','Young post','Elder pre','Elder post'};
xtickangle(45);
ylabel('Normalised variability');
ylim([0,0.3]);
yl = ylim;
for x3 = 1:size(data,2)
    p = nanmean(data(:,x3));
    y3 = yl(1)+yl(2)*0.1;%-0.1*(yl(2)-yl(1));
    txt3 = [sprintf('%0.2f',p)];
    text(x3,y3,txt3,'FontSize',9,'HorizontalAlignment','center','Rotation',90);
end

a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Arial','fontsize',8);


%% PDR evoked by a harmonic tone
pupilunit = '%'; divideByResting = 0;

load([path_in '\loc_Elder' '[' pupilunit ']'],'pupildiameter1','pupildiameter_std1','pupildiameter2','pupildiameter_std2','PDR','timeaxis','sublist');
p_e_pre = pupildiameter1(:,1)';
p_e_post = pupildiameter2(:,1)';
s_e_pre = pupildiameter_std1(:,1)';
s_e_post= pupildiameter_std2(:,1)';
pdr_e = PDR;

load([path_in '\loc_Young' '[' pupilunit ']'],'pupildiameter1','pupildiameter_std1','pupildiameter2','pupildiameter_std2','PDR','timeaxis','sublist');
p_y_pre = pupildiameter1(:,1)';
p_y_post = pupildiameter2(:,1)';
s_y_pre = pupildiameter_std1(:,1)';
s_y_post= pupildiameter_std2(:,1)';
pdr_y = PDR;

%% PDR
subplot(5,1,2); hold on;

PDRC = {};
test = [];
for g = 1:2
    if g==1
        load([path_in '\loc_Young' '[' pupilunit ']'],'pupildiameter1','pupildiameter_std1','pupildiameter2','pupildiameter_std2','PDR','timeaxis','sublist');
    elseif g==2
        load([path_in '\loc_Elder' '[' pupilunit ']'],'pupildiameter1','pupildiameter_std1','pupildiameter2','pupildiameter_std2','PDR','timeaxis','sublist');
    end
    
    if divideByResting
        p_pre = pupildiameter1(:,1);
        pdr = squeeze(nanmean(PDR,2));
        pdr = pdr./p_pre;
    else
        pdr = squeeze(nanmean(PDR,2));
    end
    
    PDRC{g} = pdr;
    
    %% Smooth
    for s=1:size(pdr,1)
        pdr(s,:) = smooth(pdr(s,:),150,'hann'); %Smooth pupil data
    end
    a = nanmean(pdr,1);
    b = nanstd(pdr)/sqrt(numel(sublist));
    
    curve1 = a'+b';
    curve2 = flipud(a'-b');
    X = [timeaxis'; flipud(timeaxis')];
    Y = [curve1; curve2];
    figfill = fill(X',Y',colourmap_group(g,:),'edgecolor','none','facealpha',0.2);
    plot(timeaxis,a,'LineWidth',1.5,'Color',colourmap_group(g,:));
    
    
    tw_analysis = [0,3];
    tx = [find(timeaxis>=tw_analysis(1)), find(timeaxis<= tw_analysis(2))];
    for s = 1:size(pdr,1)
        test(g,s) = nanmean(pdr(s,tx));
    end
end

switch pupilunit
    case 'mm'
        ylim([-0.05,0.3]);
        if divideByResting
            ylim([-0.01,0.04]);
        end
    case '%'
        ylim([0.99,1.04]);
end

ab = ylim;
aa= 0.025*(max(ab)-min(ab));
dist = 0.05*(max(ab)-min(ab));
Ystat = ab(1)+aa;
Ystat = Ystat+ dist;
cond1 = squeeze(PDRC{1});
cond2 = squeeze(PDRC{2});
dataB = bootstrap(cond1)-bootstrap(cond2);
s_boot = findSigDiff(dataB, 0.05);

plot(timeaxis,Ystat*abs(s_boot),'LineWidth',2.5,'Color','k');
%     text(0.2,ab(2),['min cluster = ',num2str(max(statlimit)*(abs(timeaxis(1))+abs(timeaxis(end)))/size(timeaxis,2)*100,'%0.2f'),' ms'],'Color','k','FontSize',12,'VerticalAlignment','middle','HorizontalAlignment','left');

hold off;
ylabel(['Pupil diameter relative to baseline' '[' pupilunit ']']);

%% PDR Derivative
subplot(5,1,5);
hold on;
D = {};
for g = 1:2
    if g==1
        load([path_in '\loc_Young' '[' pupilunit ']'],'pupildiameter1','pupildiameter_std1','pupildiameter2','pupildiameter_std2','PDR','timeaxis','sublist');
    elseif g==2
        load([path_in '\loc_Elder' '[' pupilunit ']'],'pupildiameter1','pupildiameter_std1','pupildiameter2','pupildiameter_std2','PDR','timeaxis','sublist');
    end
    p_pre = pupildiameter1(:,1);
    d = [];
    fs = 1/(timeaxis(2)-timeaxis(1));
    for s = 1:size(PDR,1) %sublist
        for i = 1:size(PDR,2) %number of trial
            pn = squeeze(PDR(s,i,:));
            pn = smooth(pn,150,'hann');
            if divideByResting
                pn = pn/p_pre(s);
            end
            d(s,i,:) = [0,diff(pn)*fs];
        end
    end
    
    derv = squeeze(nanmean(d,2));
    D{g} = derv;
    
    a = nanmean(derv,1);
    b = nanstd(derv)/sqrt(numel(sublist));
    
    curve1 = a'+b';
    curve2 = flipud(a'-b');
    X = [timeaxis'; flipud(timeaxis')];
    Y = [curve1; curve2];
    figfill = fill(X',Y',colourmap_group(g,:),'edgecolor','none','facealpha',0.2);
    plot(timeaxis,a,'LineWidth',1.5,'Color',colourmap_group(g,:));
    
end
switch pupilunit
    case 'mm'
        ylim([-0.0004,0.00075]);
        if divideByResting
            ylim([-0.000075,0.00015]);
        end
    case '%'
        ylim([-0.045,0.105]);
end

ab = ylim;
aa = 0.025*(max(ab)-min(ab));
dist = 0.05*(max(ab)-min(ab));
Ystat = ab(1)+aa;

cond1 = squeeze(D{1});
cond2 = squeeze(D{2});
dataB = bootstrap(cond1)-bootstrap(cond2);
s_boot = findSigDiff(dataB, 0.05);
plot(timeaxis,Ystat*abs(s_boot),'LineWidth',2.5,'Color','k');

hold off;

xlabel('Time from onset [s]');
ylabel(['Velocity in pupil diameter ' '[' pupilunit ']']);
% title(['PDR derivative to a harmonic tone'],'FontSize',14);

%% Compute cv on non-baseline corrected data
pupilunit = '%';
bc = 0;
filename = sprintf('%s/result_ht_%s_bc%d',path_in,pupilunit,bc);

%% Extract PDR for each group
load(filename,'PDR','cv_within','timeaxis');


%% BETWEEN subject variability: coefficient of variation
groups = {'Young','Elder'};

subplot(5,1,3);hold on;

CV_BETWEEN = [];
for group = 1:numel(groups)
    pdr = PDR{group};
    indmean = squeeze(nanmean(pdr,2));
    M{group} = indmean;
    
    indmean = M{group};
    cv = std(indmean)./(mean(indmean));
    
    CV_BETWEEN(group,:) = cv;
    
    a = cv;
    plot(timeaxis,a,'LineWidth',1.5,'Color',colourmap_group(group,:));
end
hold off;
ylabel('Between subject variability');
ylim([0.05,0.2]);

%% WITHIN SUBJECT VARIABILITY
subplot(5,1,4); hold on;

for group = 1:numel(groups)
    indmean = [];
    indstd = [];
    
    cv = cv_within{group};
    
    % PLOT
    a = nanmean(cv);
    plot(timeaxis,a,'LineWidth',1.5,'Color',colourmap_group(group,:));
end

% PLOT error bar
for group = 1:numel(groups)
    cv = cv_within{group};
    a = nanmean(cv)';
    b = (nanstd(cv)/sqrt(size(cv,1)))';
    
    curve1 = a+b;
    curve2 = flipud(a-b);
    
    X = [timeaxis'; flipud(timeaxis')];
    Y = [curve1; curve2];
    figfill = fill(X,Y,colourmap_group(group,:),'edgecolor','none','facealpha',0.2);
end
set(gca,'children',flipud(get(gca,'children'))); % Send shaded areas in background
ylabel('Within subject variability');

ylim([0.01 0.12]);

% PLOT BOOTSTRAP STATS
ab = ylim;
aa= 0.025*(max(ab)-min(ab));
dist = 0.09*(max(ab)-min(ab));

Ystat = ab(1)+aa;
cond1 = squeeze(cv_within{1});
cond2 = squeeze(cv_within{2});
dataB = bootstrap(cond1)-bootstrap(cond2);
s_boot = findSigDiff(dataB, 0.05);
plot(timeaxis,Ystat*abs(s_boot),'LineWidth',5,'Color','k');
hold off;

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 12]);

%% Statistical test in paper
% the average pupil diameter over the first 3 seconds following tone onset was significantly above floor as confirmed by a one-sample t-test in each group (young: t(17) = 5.07, p<.001; older: t(17) = 5.13, p<.001).
[h,p,ci,stats] = ttest2(test(1,:),test(2,:));

[h,p1,ci,stats1] = ttest(test(1,:),0);
[h,p2,ci,stats2] = ttest(test(2,:),0);
fprintf('t-test young: t(%d) = %.2f, p=%d; older: t(%d) = %.2f, p=%d\n',stats1.df,stats1.tstat,p1,stats2.df,stats2.tstat,p2);
