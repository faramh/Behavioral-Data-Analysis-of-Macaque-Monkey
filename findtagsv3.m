% find stim on and result tags
clc; clear all; close all;
% load data
filename = sprintf('Trial_ID_11038_Search.mat');
load(filename);
% make events to string
events = string(edf_data.Events.Messages.info(1,:));
% load Ref mat to check tags
load REF.mat
% use SearchTaskParam_v03 function
SearchEvents = SearchTaskParam_v03(edf_data,Ref);
%% build MyStruct with stimon to result epochs for each trial
% find indexes of each delimiter for trials
temp_delimiter = find(cellfun(@(x) strcmp(x,'========== STATS ========='),events));
% get number of trilas
MyStruct2.NumTrial= length(temp_delimiter) + 1;

%% indexes of stim on and result
% find stim on of each trial 
for i=1:MyStruct2.NumTrial
    temp_tag_idx=find(cellfun(@(y) isempty(y),regexp([cellfun(@(x) (x), events(SearchEvents.TrialScopeIdxs(i):SearchEvents.TrialScopeIdxs(i+1)),'UniformOutput',false)], 'TAG: [A-Z]{3}'))==0);
     stim_tag_meseg = {}; 
     stim_tag_asc = [];
    if (~isempty(temp_tag_idx))
        MyStruct2.Trial{i}.indxStimOn = temp_tag_idx + SearchEvents.TrialScopeIdxs(i)-1;
        stim_tag_meseg = events(SearchEvents.TrialScopeIdxs(i) + temp_tag_idx(1) - 1);
        stim_tag_asc = char(cellfun(@(x) x(7:end-1), stim_tag_meseg, 'UniformOutput', false));
        MyStruct2.Trial{i}.tagStimOn = (stim_tag_asc(1) - 65) * (2 ^ 10) + (stim_tag_asc(2) - 65) * (2 ^ 5) + (stim_tag_asc(3) - 65); 
        if Ref(find(MyStruct2.Trial{i}.tagStimOn==Ref(:,1),1),4)==2
            MyStruct2.Trial{i}.istagStimOn = 1;
        end
    end
end
% find Reuslt tag in each trial scope
for i=1:MyStruct2.NumTrial
    MyStruct2.Trial{i}.indxResult = NaN;
    temp_tag_idx = find(cellfun(@(y) isempty(y), regexp([cellfun(@(x) (x), events(SearchEvents.TrialScopeIdxs(i):SearchEvents.TrialScopeIdxs(i+1)),'UniformOutput',false)],'Result Enter'))==0);
    if ~isempty(temp_tag_idx)
        MyStruct2.Trial{i}.indxResult = temp_tag_idx + SearchEvents.TrialScopeIdxs(i)-1;
    end
end
%%
% get eye segmenter result for each trial with stim on index
if length(edf_data.Events.Eblink.eye{1}) == 4, eyeNum = 1; else eyeNum = 2; end

gx = edf_data.Samples.gx(:,eyeNum).'; gy = edf_data.Samples.gy(:,eyeNum).'; % Prepare Input
gx(isnan(gx)) = -10000; gy(isnan(gy)) = -10000; %Drop out nans
gaze{1,1} = gx; gaze{1,2} = gy; %Put into cell
false_brks_ths = 0.025; minisac_ths =0.1;
Segment=Eyepos_Segmenter02({gx, gy},{0,[],false_brks_ths,minisac_ths,[],1},'AutoPromAdjust',0);
% add time to the Segmat matrix (column 8 and 9)
Segment(:, 8) = edf_data.Samples.time(Segment(:,1))- edf_data.Samples.time(1);
Segment(:, 9) = edf_data.Samples.time(Segment(:,2))- edf_data.Samples.time(1);
%% separte esegments for each trial with stim on
for i=1:MyStruct2.NumTrial
    if isfield(MyStruct2.Trial{i},'indxStimOn')
        trl_ts = edf_data.Events.Messages.time(MyStruct2.Trial{i}.indxStimOn) - edf_data.Samples.time(1);
        trl_te = edf_data.Events.Messages.time(MyStruct2.Trial{i}.indxResult) - edf_data.Samples.time(1);
        MyStruct2.Trial{i}.Segment = Segment(Segment(:,8)>=trl_ts & Segment(:,9)<=trl_te,:);
    end
end
%% plot fractal positions and eye fixed points for all trials 

k=1;
for i=1:MyStruct2.NumTrial
    if (isfield(MyStruct2.Trial{i},'indxStimOn') & SearchEvents.TrialsProperties{i}.TrialValidity==1)
        for j=1:SearchEvents.TrialsProperties{i}.DisplaySize
            fx(k)=SearchEvents.TrialsProperties{i}.Region(j).x;
            fy(k)=SearchEvents.TrialsProperties{i}.Region(j).y;
            k=k+1;
        end
    end
end

% gx and gy for fixated points
k=1;
for i=1:MyStruct2.NumTrial
    if (isfield(MyStruct2.Trial{i},'indxStimOn') & SearchEvents.TrialsProperties{i}.TrialValidity==1 & SearchEvents.TrialsProperties{i}.Accuracy==1)
        for j=1:size(MyStruct2.Trial{i}.Segment,1)
            if (MyStruct2.Trial{i}.Segment(j,5)==2 & MyStruct2.Trial{i}.Segment(j,3)~=-10000)
                gx_fix(k)=MyStruct2.Trial{i}.Segment(j,3);
                gy_fix(k)=MyStruct2.Trial{i}.Segment(j,4);
                k=k+1;
            end
        end
    end
end
f1=figure(1);
plot(fx,fy,"*",'Color','r','MarkerSize',9,'LineWidth',1)
hold on
plot(gx_fix,gy_fix,"x",'Color','b','MarkerSize',9,'Linewidth',1)
hold 
xlim([-500 1500]);
ylim([-500 1500]);
axis equal
grid on
title('fractal and gaze positions before tranformation')
savefilename=sprintf('beforetransition_v2');
print(f1,savefilename,'-dpng','-r300');
%% find transformatiuon parameter
% l=1;
% 
% data(:,1)=gx_fix;
% data(:,2)=gy_fix;
% 
% % Specify the number of clusters you want to find
% numClusters =1;
% % Perform K-means clustering
% [idx, C] = kmeans(data, numClusters);
% % Plot the data points with different colors for each cluster
% figure;
% gscatter(data(:,1), data(:,2), idx);
% title('K-means Clustering');
% % Plot the cluster centers
% hold on;
% plot(C(:,1)+15, C(:,2)-20, 'kx', 'MarkerSize', 15, 'LineWidth', 3);
% xlim([0 900]);
% ylim([0 900]);
% axis equal
% %finding trnsofmation value
% %fit a circle to fractal positions
% [xc, yc, R]=lsqcircle(fx', fy');
% center=[xc, yc]; % Center coordinates (x, y)
% radius=R;      % Radius of the circle
% dx=round(C(:,1)+15-xc,0);
% dy=round(C(:,2)-10-yc,0);
%%
dx=210;
dy=140;
f3=figure(3);
plot(fx,fy,"*",'Color','r','MarkerSize',9,'LineWidth',1)
hold on
plot(gx_fix-dx,gy_fix-dy,"x",'Color','b','MarkerSize',9,'Linewidth',1)

xlim([0 1500]);
ylim([0 1500]);
axis equal
title('fractal and gaze positions after tranformation')
savefilename=sprintf('afterv2');
%% plot eye trajectory for correct TP trials
% find TP correct index 
indxTPCorr=[];
for i=1:MyStruct2.NumTrial
    if (isfield(MyStruct2.Trial{i},'indxStimOn') & SearchEvents.TrialsProperties{i}.TargetPresence==1 &...
       SearchEvents.TrialsProperties{i}.Accuracy==1)
        indxTPCorr=[indxTPCorr, i];
    end
end

% plot five TP correct trials
f4=figure(4);
%subplot(5,1,1)
%fractals
n=6;
fractal_win= 144; % size of fractal window
for i=1:SearchEvents.TrialsProperties{indxTPCorr(n)}.DisplaySize
    x1=SearchEvents.TrialsProperties{indxTPCorr(n)}.Region(i).x;
    y1=SearchEvents.TrialsProperties{indxTPCorr(n)}.Region(i).y;
    rectangle('Position', [x1, y1, fractal_win, fractal_win], 'EdgeColor', 'r');
    axis equal;
    hold on
end
% transformed gx and gy before first saccade
%temp_idx=find(MyStruct.Trial{indxTPCorr(n)}.Segment(:,5)==2,1);

% plot for the whole epoch of trial

% find time of sample form SearchParam output
t_sac_ini = edf_data.Samples.time(MyStruct2.Trial{indxTPCorr(n)}.idx_ts + MyStruct2.Trial{indxTPCorr(n)}.Segment(1,1));
t_sac_lan = edf_data.Samples.time(MyStruct2.Trial{indxTPCorr(n)}.idx_ts + MyStruct2.Trial{indxTPCorr(n)}.Segment(end,2));
idx_temp_ini=find(edf_data.Samples.time(:,1)==t_sac_ini,1);
idx_temp_lan=find(edf_data.Samples.time(:,1)==t_sac_lan,1);
gx_trans=edf_data.Samples.gx(idx_temp_ini:idx_temp_lan)- dx;
gy_trans=edf_data.Samples.gy(idx_temp_ini:idx_temp_lan)- dy;
colors = jet(size(gx_trans, 2));

data_x = gx_trans(1:2:end);
data_y = gy_trans(1:2:end);
scatter(data_x, data_y, [], colors(1:length(data_x),:,:), 'filled');
%plot(gx_trans,gy_trans,"x",'Color','b','MarkerSize',9,'Linewidth',1)



