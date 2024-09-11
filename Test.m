% This code is arranged to start working with behvaioral data 
%% loading the data
clc; clear; close all;
% Add data filder with subfolders
addpath(genpath('Data'));
addpath(genpath('FC Task'));



%% Data Extraction

load("Data/Trial_ID_9802_FC.mat")
load("REF.mat");

% Using FCTaskParam_v03 function
% We Extraxt parameters of FC session and arrange them for trials 
FCEvents = FCEVENTExtractor(edf_data,Ref);

    
Trials_Table=FCEvents.TrialsProperties; %Ectraxt trial Props
TrialNumber=43; % A number between 1 to NumTrials
CurrentTrial=Trials_Table(1,TrialNumber);
NextTrial=Trials_Table(1,TrialNumber+1); % To access ITI
TrialType=CurrentTrial.type;
%% Extracting Main part of Trial

% if TrialType=='Force'
Event_Data=CurrentTrial.data;
% Example Event_Data (you can replace this with your actual data)
% Event_Data = { 'event_tag1', time1; 'event_tag2', time2; ... }


Tag_Enter='FractOverlap Enter';
% Find indices for the event "looked into fixation window"
idx_looked_into_fixation = find(contains(Event_Data(:,1), Tag_Enter));
Tag_Exit='FractFix Exit';
% Find indices for the event "FractFix Exit"
idx_fractfix_exit = find(contains(Event_Data(:,1), Tag_Exit));

% Extract the corresponding times
TrialEnter = Event_Data{idx_looked_into_fixation, 2};
TrialExit = Event_Data{idx_fractfix_exit, 2};


%% Plot the gaze in the main part of trial


% Define parameters
dx = 160;
dy = 80;
fractal_win = 144;
fixation_x = 460;
fixation_y = 390;
fixation_w = 40;
fixation_h = 40;



%% PLOT after trial ITI duration

NextEvent_Data=NextTrial.data;
% Example Event_Data (you can replace this with your actual data)
% Event_Data = { 'event_tag1', time1; 'event_tag2', time2; ... }


Tag_Enter_ITI='ITI Enter';
% Find indices for the event "looked into fixation window"
idx_Entered_ITI = find(contains(NextEvent_Data(:,1), Tag_Enter_ITI));
Tag_Exit_ITI='ITI Exit';
% Find indices for the event "FractFix Exit"
idx_exited_ITI = find(contains(NextEvent_Data(:,1), Tag_Exit_ITI));

% Extract the corresponding times
TrialEnterITI = NextEvent_Data{idx_Entered_ITI, 2};
TrialExitITI = NextEvent_Data{idx_exited_ITI, 2};




%% 
% Create a new figure for both plots
figure;

% First subplot for the gaze trace
subplot(1,2,1);  % 2 rows, 1 column, plot 1
plotGazeTrace(CurrentTrial, TrialEnter, TrialExit, dx, dy, fractal_win, fixation_x, fixation_y, fixation_w, fixation_h);
title('Gaze Trace During Trial', 'FontSize', 14, 'FontWeight', 'bold');

% Second subplot for ITI duration
subplot(1,2,2);  % 2 rows, 1 column, plot 2
PlotITIGazeTrace(NextTrial, CurrentTrial, TrialEnterITI, TrialExitITI, dx, dy, fractal_win, fixation_x, fixation_y, fixation_w, fixation_h);
title('Gaze Trace During ITI', 'FontSize', 14, 'FontWeight', 'bold');






