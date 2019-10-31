clear;close all;

% experimentname = 'Expt1 - Young'; % Uncomment to plot figure 2
experimentname = 'Expt2 - Old'; % Uncomment to plot figure 6

path_root = 'E:\Sijia\EALOPs\Submission\Upload\Processed data_Main';
path_in = [path_root '\' experimentname];

blocklist = 1:6;
trial_number = 12;
condlist = {'Easy','Medium','Hard'};%{'1 stream','2 streams','3 streams'};

%% Retrieve subject informtion from excel file
[~,~,expMat] = xlsread([path_in '\subject_Info.xlsx']);
FsMatrix = expMat(:,4); %sampling frequency, you can set it as 1000Hz
sublist = expMat(:,1); %cell array for subject names
smpfreq = FsMatrix{1};
define_colourmap;

filename = [path_in, '\behav_result.mat'];
load(filename);

filename = [path_in, '\sublist_group.mat'];
load(filename,'sublist_good','sublist_poor','sublist_exclude','trial2keep_valid');
sublist_include = [sublist_good;sublist_poor];

%% All performers
figure(1);clf;
subplot(3,1,1);
data = avg_pHit*100;
plotbar(data,0);
ax=gca;
ax.XTick = 1:1:size(data,2);
ax.XTickLabel = condlist;
%     xtickangle(45)
% set(gca, 'YTick', 0:20:100);
ylabel('Hit rate [%]');
%     title('Pupil diameter at resting');
set(gca,'FontSize',12);
yl = [40 100]; %[mm]
ylim(yl);
for x3 = 1:size(data,2)
    p = nanmean(data(:,x3));
    y3 = yl(1)+yl(2)*0.1;%-0.1*(yl(2)-yl(1));
    txt3 = [sprintf('%0.2f',p)];
    text(x3,y3,txt3,'FontSize',9,'HorizontalAlignment','center','Rotation',90);
end

subplot(3,1,2);
data = sum_nFa;
plotbar(data,0);
ax=gca;
ax.XTick = 1:1:size(data,2);
ax.XTickLabel = condlist;
%     xtickangle(45)
% set(gca, 'YTick', 0:20:100);
ylabel('#FA');
%     title('Pupil diameter at resting');
set(gca,'FontSize',12);
yl = [0 20]; %[mm]
ylim(yl);
for x3 = 1:size(data,2)
    p = nanmean(data(:,x3));
    y3 = yl(1)+yl(2)*0.1;%-0.1*(yl(2)-yl(1));
    txt3 = [sprintf('%0.2f',p)];
    text(x3,y3,txt3,'FontSize',9,'HorizontalAlignment','center','Rotation',90);
end

subplot(3,1,3);
data = num_badtrials;
plotbar(data,0);
ax=gca;
ax.XTick = 1:1:size(data,2);
ax.XTickLabel = condlist;
%     xtickangle(45)
% set(gca, 'YTick', 0:20:100);
ylabel('#Bad trials');
%     title('Pupil diameter at resting');
set(gca,'FontSize',12);
yl = [0 24]; %[mm]
ylim(yl);
for x3 = 1:size(data,2)
    p = nanmean(data(:,x3));
    y3 = yl(1)+yl(2)*0.1;%-0.1*(yl(2)-yl(1));
    txt3 = [sprintf('%0.2f',p)];
    text(x3,y3,txt3,'FontSize',9,'HorizontalAlignment','center','Rotation',90);
end

suptitle(['All participants (N=' num2str(size(data,1)) ')']);
