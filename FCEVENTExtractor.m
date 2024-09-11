function FCEvents = FCEVENTExtractor(edf_data,Ref)

% -------------------------------------------------------------------------
%  Made By Mohamad Hosein Faramarzi 2024
%-------------------------------------------------------------------------- 

MassageTime = edf_data.Events.Messages.time';
MassageInfo = edf_data.Events.Messages.info';
MassageInfoStr = convertCharsToStrings(MassageInfo);

% Convert MassageTime to a cell array
MassageTimeCell = num2cell(MassageTime);

% Concatenate MassageTimeCell and MassageInfo
FullMessage2 = [MassageInfoStr,MassageTimeCell];

FullMessage = [MassageInfo,MassageTimeCell];

Margin_Insufficient_Trial_Number=15;
disp('Evaluating FC task')
% Initialize parameters
FCEvents = struct();

FCEvents.FullMessageInfo=FullMessage2;



% Define the tags that mark the start of trials
tagForce = 'Current Experiment : Force';
tagChoice = 'Current Experiment : Choice';

% Find the indices of the start tags
startIndicesForce = find(strcmp(FullMessage(:,1), tagForce));
startIndicesChoice = find(strcmp(FullMessage(:,1), tagChoice));

% Combine the indices and sort them to get the trial start points in order
startIndices = sort([startIndicesForce; startIndicesChoice]);

% Define the tag for ITI
itiTag = '-iti ';
setTag= '-set ';

% Initialize a struct to store trials
Trials = struct();

% Loop through the indices and extract trials
for i = 1:length(startIndices)
    
    % Get the starting index for this trial
    startIdx = startIndices(i);
 
    % Determine the trial type based on the tag
    startTag = FullMessage{startIdx, 1};
    if contains(startTag, tagForce)
        Trials(i).type = 'Force';
    elseif contains(startTag, tagChoice)
        Trials(i).type = 'Choice';
    else
        Trials(i).type = 'Unknown'; % Fallback if type is not identified
    end

    % Determine the end index of the trial
    if i < length(startIndices)
        endIdx = startIndices(i+1) - 1; % Go up to the next start tag
    else
        endIdx = size(FullMessage, 1); % For the last trial, go to the end
    end
    
    % Extract the data for this trial
    Trials(i).data = FullMessage(startIdx:endIdx, :);
    Trials(i).FullMessage=FullMessage2(startIdx:endIdx, :);
    
    % Look for the ITI tag within this trial
    itiIdx = find(contains(FullMessage(startIdx:endIdx, 1), itiTag));
    
    if ~isempty(itiIdx)
        % If ITI tag is found, extract the corresponding number using regexp
        itiString = FullMessage{startIdx + itiIdx(1) - 1, 1}; % Get the string with -iti
        
        % Use regexp to extract the number following '-iti'
        itiNumber = regexp(itiString, '-iti\s*(\d+)', 'tokens');
        
        if ~isempty(itiNumber)
            % Convert the extracted number (which is a string) to a number
            Trials(i).iti = str2double(itiNumber{1}{1});
        else
            % If no number found, set ITI to NaN
            Trials(i).iti = NaN;
        end
    else
        % If no ITI tag is found, set ITI to NaN
        Trials(i).iti = NaN;
    end




    setIdx = find(contains(FullMessage(startIdx:endIdx, 1), setTag));


    if ~isempty(setIdx)
        % If SET tag is found, extract the corresponding number using regexp
        setString = FullMessage{startIdx + setIdx(1) - 1, 1}; % Get the string with -set
        
        % Use regexp to extract the number following '-set'
        setNumber = regexp(setString, '-set\s*(\d+)', 'tokens');
        
        if ~isempty(setNumber)
            % Convert the extracted number (which is a string) to a number
            Trials(i).set = str2double(setNumber{1}{1});
        else
            % If no number found, set SET to NaN
            Trials(i).set = NaN;
        end
    else
        % If no SET tag is found, set SET to NaN
        Trials(i).set = NaN;
    end


end
% Store the number of trials in a variable
numTrials = length(startIndices);

% Store the data to main struct
FCEvents.NumTrials = length(startIndices);


if (FCEvents.NumTrials < Margin_Insufficient_Trial_Number)
    FCEvents.Invalidated = 1;
    FCEvents.Ref = [];
    warning(['Number of trials was not enough for task #'])
    return;
else
    FCEvents.Invalidated = 0;
end
%% 
% Eye properties:
% left or right eye:
if length(edf_data.Events.Eblink.eye{1}) == 4, eyeNum = 1; else eyeNum = 2; end

gx = edf_data.Samples.gx(:,eyeNum).'; gy = edf_data.Samples.gy(:,eyeNum).'; % Prepare Input
gx(isnan(gx)) = -10000; gy(isnan(gy)) = -10000; %Drop out nans
gaze{1,1} = gx; gaze{1,2} = gy; %Put into cell


% Segment
false_brks_ths = 0.025; minisac_ths =0.1;

[Segmat,Eye_vel] = Eyepos_Segmenter02(gaze,{0,[],false_brks_ths,minisac_ths,[],1}, ...
    'AutoPromAdjust',0);

%Drop out Nans
Segmat(...
    Segmat(:,3) == -10000 & ... %If it landed on -10000
    Segmat(:,5) == 0 ... % And the type was stationary
    ,:) = [];

% add time to the Segmat matrix (column 8 and 9)
Segmat(:, 8) = edf_data.Samples.time(Segmat(:,1))- edf_data.Samples.time(1);
Segmat(:, 9) = edf_data.Samples.time(Segmat(:,2))- edf_data.Samples.time(1);
FCEvents.Segmat = Segmat;



% Initialize gaze data fields in the Trials struct
for i = 1:length(Trials)
    Trials(i).gx = [];
    Trials(i).gy = [];
    Trials(i).time = [];
end

% Loop through each trial and extract gaze data
for i = 1:length(Trials)
    % Extract the start and end times of the trial
    startTime = Trials(i).data{1, 2};  % First row, second column (time)
    endTime = Trials(i).data{end, 2};  % Last row, second column (time)
    
    % Find the index range in edf_data.Samples.time that corresponds to the trial
    startIdx = find(edf_data.Samples.time >= startTime, 1, 'first');
    endIdx = find(edf_data.Samples.time <= endTime, 1, 'last');
    
    % Extract gaze positions (gx and gy) for the trial
    if ~isempty(startIdx) && ~isempty(endIdx)
            % gaze epochs


        Trials(i).gx = edf_data.Samples.gx(startIdx:endIdx,2);
        Trials(i).gy = edf_data.Samples.gy(startIdx:endIdx,2);

        Trials(i).time = edf_data.Samples.time(startIdx:endIdx);
        Trials(i).gxvel = edf_data.Samples.gxvel(startIdx:endIdx,2);
        Trials(i).gyvel = edf_data.Samples.gyvel(startIdx:endIdx,2);
        %Segmenter output Abbaszade version
        Trials(i).segmat=Segmat(Segmat(:,8)>=startIdx & Segmat(:,9)<=endIdx,:);

        % drop nans
        Trials(i).gx(isnan(Trials(i).gx)) = -10000; Trials(i).gy(isnan(Trials(i).gy)) = -10000;
        if ~isempty(find(Trials(i).gx ~= -10000))
        % eye segmenter function
            false_brks_ths = 0.025; minisac_ths =0.1; 
            concat={Trials(i).gx,Trials(i).gy};
            Gx=Trials(i).gx.';
            Gy=Trials(i).gy.';
            Gaze{1,1}=Gx; Gaze{1,2}=Gy;
            Trials(i).NewSegment = Eyepos_Segmenter02(Gaze,{0,[],false_brks_ths,minisac_ths,[],1},'AutoPromAdjust',0);
        end


        Trials(i).NewSegment(:, 8) = Trials(i).time(Trials(i).NewSegment(:,1)); % - edf_data.Samples.time(1);
        Trials(i).NewSegment(:, 9) = Trials(i).time(Trials(i).NewSegment(:,2)); % - edf_data.Samples.time(1);
    else
        % Handle the case where no matching time points were found (optional)
        Trials(i).gx = NaN;
        Trials(i).gy = NaN;
        Trials(i).time = NaN; 
        Trials(i).gxvel = NaN;
        Trials(i).gyvel = NaN;   
        Trials(i).segmat = NaN;

    end
end

%% 

for i = 1:length(Trials)
    % Initialize default values for trial type, validity, result, and correctness
    Trials(i).TrialType = -1; % -1 indicates unknown trial type
    Trials(i).TrialValidity = 0; % 0 means invalid
    Trials(i).TrialCorrectness = 0; % 0 means incorrect
    Trials(i).Reward = -1; % -1 indicates unknown result
    
    % Extract the trial event data for this trial
    trialEvents = Trials(i).data(:, 1); % First column contains event labels
    
    % Determine the trial type
    if any(contains(trialEvents, 'Current Experiment : Force'))
        Trials(i).TrialType = 0; % 0 for 'Force'
    elseif any(contains(trialEvents, 'Current Experiment : Choice'))
        Trials(i).TrialType = 1; % 1 for 'Choice'
    end
    
    % Validate the trial based on certain criteria
    if Trials(i).TrialType == 0 % Force trial
        if any(contains(trialEvents, 'Sound Success Played'))
            Trials(i).TrialValidity = 1; % Valid trial
            
            % Check the trial result for Force trials
            if any(contains(trialEvents, 'Good Fractal Reward Delivered'))
                Trials(i).Reward = 'Good'; % Correct result
                Trials(i).TrialCorrectness = 1; % Correct trial
            elseif any(contains(trialEvents, 'Bad Fractal Reward Delivered'))
                Trials(i).Reward = 'Bad'; % Incorrect result
                Trials(i).TrialCorrectness = 0; % Incorrect trial
            else
                Trials(i).Reward = -1; % Unknown result
                Trials(i).TrialCorrectness = 0; % Unknown correctness
                warning('Suspicious Force trial found in Trial %d', i);
            end
        end
        
    elseif Trials(i).TrialType == 1 % Choice trial
        if any(contains(trialEvents, 'Sound Success Played'))
            Trials(i).TrialValidity = 1; % Valid trial
            
            % Check the trial result for Choice trials
            if any(contains(trialEvents, 'Good Fractal Reward Delivered'))
                Trials(i).Reward = 'Good'; % Correct result
                Trials(i).TrialCorrectness = 1; % Correct trial
            elseif any(contains(trialEvents, 'Bad Fractal Reward Delivered'))
                Trials(i).Reward = 'Bad'; % Incorrect result
                Trials(i).TrialCorrectness = 0; % Incorrect trial
            else
                Trials(i).Reward = -1; % Unknown result
                Trials(i).TrialCorrectness = 0; % Unknown correctness
                warning('Suspicious Choice trial found in Trial %d', i);
            end
        end
    else
        % Handle the case where the trial type is unknown
        warning('Unknown trial type found in Trial %d', i);
    end
end
%% 

for i = 1:length(Trials)
    % Initialize default fields for fractal information
    Trials(i).FractalRegions = []; % Will store X and Y for Force trials or N and Z for Choice trials
    
    % Extract the trial event data for this trial
    trialEvents = Trials(i).data(:, 1); % First column contains event labels
    
    % Check trial type
    if Trials(i).TrialType == 0 % Force trial
        % Find the event with 'Fractal X is set for region Y'
        fractal_idx = find(~cellfun('isempty', regexp(trialEvents, 'Fractal \d+ is set for region \d+ with x,y,w,h =')));
        
        if ~isempty(fractal_idx)
            % Extract the event string
            fractal_event = trialEvents{fractal_idx(1)};
            
            % Extract X and Y using regexp
            fractal_info = regexp(fractal_event, 'Fractal (\d+) is set for region (\d+)', 'tokens');
            if ~isempty(fractal_info)
                X = str2double(fractal_info{1}{1});
                Y = str2double(fractal_info{1}{2});
                % Store X and Y in the Trials struct 'Fractal X is set for region Y'
                Trials(i).FractalRegions = Y;
            end
        end
        
    elseif Trials(i).TrialType == 1 % Choice trial
        % Find the event with 'Good Fractal M is set for region N'
        good_fractal_idx = find(~cellfun('isempty', regexp(trialEvents, 'Good Fractal \d+ is set for region \d+ with x,y,w,h =')));
        
        % Find the event with 'Bad Fractal H is set for region Z'
        bad_fractal_idx = find(~cellfun('isempty', regexp(trialEvents, 'Bad Fractal \d+ is set for region \d+ with x,y,w,h =')));
        
        if ~isempty(good_fractal_idx) && ~isempty(bad_fractal_idx)
            % Extract the good fractal event string
            good_fractal_event = trialEvents{good_fractal_idx(1)};
            bad_fractal_event = trialEvents{bad_fractal_idx(1)};
            
            % Extract N (region for Good Fractal) using regexp
            good_fractal_info = regexp(good_fractal_event, 'Good Fractal \d+ is set for region (\d+)', 'tokens');
            % Extract Z (region for Bad Fractal) using regexp
            bad_fractal_info = regexp(bad_fractal_event, 'Bad Fractal \d+ is set for region (\d+)', 'tokens');
            
            if ~isempty(good_fractal_info) && ~isempty(bad_fractal_info)
                N = str2double(good_fractal_info{1}{1}); % Region N for Good Fractal
                Z = str2double(bad_fractal_info{1}{1});  % Region Z for Bad Fractal
                % Store N and Z in the Trials struct
                Trials(i).FractalRegions = [N, Z];
            end
        end
    end
end

% Fractal Type extraction
for i = 1:length(Trials)
    if strcmp(Trials(i).type, 'Force') && any(Trials(i).FractalRegions < 4)
        Trials(i).FractalType = 'Good';
    elseif strcmp(Trials(i).type, 'Force') && any(Trials(i).FractalRegions >= 4)
        Trials(i).FractalType = 'Bad';
    else
        Trials(i).FractalType = 'Good,Bad';
    end
end



FCEvents.TrialsProperties=Trials;
FCEvents.Samples=edf_data.Samples;
FCEvents.REF=Ref;
end