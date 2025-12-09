%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            SHINE LAB, CU BOULDER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the main analysis script for fNIRS neuroimaging data analysis.
% 2024 - Institute of Cognitive Science, CU-Boulder - James Crum  
% 2025 - Modifications Jalynn Nicoly

%% Add NIRS Toolbox to the Matlab path
clear, clc 
addpath(genpath('C:\Program Files\MATLAB\R2023a\toolbox\nirs-toolbox-master'));

%% Initialize data directory
% Directory = 'C:\Users\Jalynn\Desktop\fNIRS_Analysis_ONR2_Germane\fNIRS';
Directory = 'C:\Users\Jalynn\Documents\GitHub\Universal_fNIRS_Analyzer\fNIRS';

% Load raw fNIRS intensities
rawDir = dir(Directory);
raw = nirs.io.loadDirectory(Directory, {'Subject'});

% Load conditions
conditions = readtable('Conditions.csv');

%% Extract, organize, & insert stimulus designs (LSL Files)
% Find all .tri files recursively
triFiles = dir(fullfile(Directory, '**', '*.tri'));
% Initialize struct to group by participant/run
trigs_by_ptp = struct();
for i = 1:numel(triFiles)
    % --- Extract participant folder (first folder under Old_fNIRS) ---
    relPath = strrep(triFiles(i).folder, [Directory filesep], '');  % remove base path
    parts = split(relPath, filesep);                                 % split remaining path
    ptpFolderName = parts{1};                                        % first folder = participant
    
    % --- Extract participant ID (P001, P002, etc.) ---
    ptpID = extractBefore(ptpFolderName, ' ');  % take part before first space
    
    % Load .tri file
    triPath = fullfile(triFiles(i).folder, triFiles(i).name);
    triData = readmatrix(triPath, 'FileType', 'text', 'Delimiter', ';');
    % Store column 3 if it exists
    trigValues = [];
    if size(triData,2) >= 3
        trigValues = triData(:,3);
    end
    % Insert into struct under participant/run
    if isfield(trigs_by_ptp, ptpID)
        error('Participant ID %s already exists in trigs_by_ptp!', ptpID);
        % trigs_by_ptp.(ptpID){end+1} = trigValues; % append to cell array
    else
        trigs_by_ptp.(ptpID) = {trigValues};       % initialize cell array
    end
end

%% Remove particiants with insufficient timestamps across vars
a = 1;
idx = [];

ptpNames = fieldnames(trigs_by_ptp);
for i = 1:numel(ptpNames)
    cellVal = trigs_by_ptp.(ptpNames{i});  % should be a 1x1 cell containing the vector
    if ~iscell(cellVal) || isempty(cellVal) || isempty(cellVal{1})
        cellSize = 0;
    else
        cellSize = size(cellVal{1},1);
    end

    if cellSize < 8
        idx(a) = i;
        a = a + 1;
    end
end

% remove collected participants (if any)
if ~isempty(idx)
    error(['Error: The following participants have fewer than 8 triggers: ' ...
            num2str(idx)]);

    % % Which ptps to remove - Those with less than 8 triggers
    % namesToRemove = ptpNames(idx);
    % 
    % % Remove ptps from trigger struct
    % trigs_by_ptp = rmfield(trigs_by_ptp, namesToRemove);
    % 
    % % Remove same ptps from raw
    % raw(idx) = [];
    % 
    % % rms same ptps from conditions dataset
    % conditions(idx,:) = [];
end

%% Convert condition codes to string labels
% Converts a table to a cell matrix
conditions = table2cell(conditions);

for i = 1:size(conditions, 1)
    for j = 2:3
        if conditions{i,j} == 1
            conditions{i,j} = 'LG';
        elseif conditions{i,j} == 2
            conditions{i,j} = 'HG';
        end
    end
end

%% Extract, organize, & insert stimulus designs (LSL Files)
% Find all .tri files recursively
triFiles = dir(fullfile(Directory, '**', '*.tri'));
% Initialize struct to group by participant/run
trigs_by_ptp = struct();

for i = 1:numel(triFiles)
    % --- Extract participant ID (P001, P002, etc.) ---
    relPath = strrep(triFiles(i).folder, [Directory filesep], '');  % remove base path
    parts = split(relPath, filesep);                                 % split remaining path
    ptpFolderName = parts{1};                                        % first folder = participant
    ptpID = extractBefore(ptpFolderName, ' ');  % take part before first space

    % Load .tri file
    triPath = fullfile(triFiles(i).folder, triFiles(i).name);

    % --- CORRECTION: Use readtable to preserve string columns ---
    % Read the file as a table, treating the delimiter correctly
    T = readtable(triPath, 'FileType', 'text', 'Delimiter', ';', 'ReadVariableNames', false);

    % Convert the table to a cell array for easy column extraction
    triData = table2cell(T);

    % Initialize variables for storage
    triTimesSec = [];
    trigValues = [];

    % Process if we have enough columns (at least timestamp and trigger code)
    if size(triData,2) >= 3
        % Extract the first column (Timestamps) as a cell array of strings
        tsColumn_Cell = triData(:, 1);

        % Convert the cell array of strings to a string array
        tsColumn = string(tsColumn_Cell);

        % Truncate microseconds to milliseconds (LSL standard)
        % Note: The strings should now be preserved, not NaN
        tsColumn = extractBefore(tsColumn, 24); % keeps yyyy-MM-ddTHH:mm:ss.SSS

        % Convert to datetime object
        triTimestamps = datetime(tsColumn, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS');

        % Calculate relative times in seconds from the first event's time
        triTimesSec = seconds(triTimestamps - triTimestamps(1));

        % Trigger codes (third column) - ensure these are converted to numbers
        trigValues = cell2mat(triData(:, 3)); 
    end

    % Store the processed data (relative times and trigger codes)
    if isfield(trigs_by_ptp, ptpID)
        error('Participant ID %s already exists in trigs_by_ptp!', ptpID);
    else
        % Store a struct containing both the time vector and trigger codes
        trigs_by_ptp.(ptpID) = struct('Times', {triTimesSec}, 'Triggers', {trigValues}); 
    end
end

%% Insert stimulus designs
%% Insert stimulus designs
nRaw = length(raw);
% Create a cell array to store the full structure from trigs_by_ptp
ptpNames = fieldnames(trigs_by_ptp);
data_trigs = cell(nRaw,1); 
for i = 1:min(nRaw,numel(ptpNames))
    % FIX: Assuming data_trigs{i} holds a struct array, index the first element
    data_trigs{i} = trigs_by_ptp.(ptpNames{i});
end

AllStims = cell(nRaw,1);
exclude = [];
a=1;

% Loop over each subject saved in raw
for i = 1:nRaw    
    tempStruct  = data_trigs{i};
    % FIX: Accessing struct fields with array indexing (1)
    if iscell(tempStruct)
         % If it somehow retrieved a cell holding the struct (as a fallback)
         tempStruct = tempStruct{1};
    end
    
    triTimesSec = tempStruct(1).Times;      % The relative time vector in seconds
    trigVec     = tempStruct(1).Triggers;   % The trigger code vector
    
    % Check for sufficient trigger points (must be >= 8: 4 starts, 4 ends)
    if isempty(trigVec) || numel(trigVec) < 8 || isempty(triTimesSec)
        exclude(a) = i; a = a+1; continue
    end
    
    % Secondary Exclusion Check (Raw NIRS Data Integrity)
    if isempty(raw(i).time) || ~isnumeric(raw(i).time)
        exclude(a) = i; a = a+1; continue
    end

    % 1. Extract Trial Onsets and Durations 
    % Pair triggers: odd = start, even = end (Indices of triTimesSec)
    starts = 1:2:length(trigVec);
    ends   = 2:2:length(trigVec);
    
    % Onsets for the 4 experimental trials 
    trial_onsets    = triTimesSec(starts);
    trial_end_times = triTimesSec(ends); 
    trial_durations = trial_end_times - trial_onsets;

    % 2. Define Condition Labels
    % C1 is conditions{i,2} (e.g., 'LG'), C2 is conditions{i,3} (e.g., 'HG')
    % Order: C1 Train (1), C1 Test (2), C2 Train (3), C2 Test (4)
    cond_name_1 = conditions{i,2};
    cond_name_2 = conditions{i,3};
    
    % Use these base names to create the explicit stimulus names
    % Order: C1 Train, C1 Test, C2 Train, C2 Test
    base_names = { cond_name_1, cond_name_1, cond_name_2, cond_name_2 };
    
    % Need 6 total events: 4 trials + 2 baselines
    newStims = cell(6,1); 
    
    % 3. Calculate Dynamic Baseline Durations
    
    % --- Baseline 1: Between Cond 1 Train (Trial 1) and Cond 1 Test (Trial 2) ---
    baseline1_onset    = trial_end_times(1);     % Starts when Trial 1 ends
    baseline1_end      = trial_onsets(2);        % Ends when Trial 2 starts
    baseline1_duration = baseline1_end - baseline1_onset;
    
    % --- Baseline 2: Between Cond 2 Train (Trial 3) and Cond 2 Test (Trial 4) ---
    baseline2_onset    = trial_end_times(3);     % Starts when Trial 3 ends
    baseline2_end      = trial_onsets(4);        % Ends when Trial 4 starts
    baseline2_duration = baseline2_end - baseline2_onset;


    % 4. Insert All 6 Stimulus Events in Order
    
    % --- Event 1: Cond 1 Train ---
    idx = 1;
    newStims{idx}.names     = [base_names{1} '_Train']; % e.g., 'LG_Train'
    newStims{idx}.onsets    = trial_onsets(1);
    newStims{idx}.durations = trial_durations(1);
    newStims{idx}.amp       = 1;
    newStims{idx}.regressor_no_interest = 0;

    % --- Event 2: Cond 1 Baseline ---
    idx = 2;
    newStims{idx}.names     = [base_names{1} '_Baseline']; % e.g., 'LG_Baseline'
    newStims{idx}.onsets    = baseline1_onset;
    newStims{idx}.durations = baseline1_duration;
    newStims{idx}.amp       = 1;
    newStims{idx}.regressor_no_interest = 0;

    % --- Event 3: Cond 1 Test ---
    idx = 3;
    newStims{idx}.names     = [base_names{2} '_Test']; % e.g., 'LG_Test'
    newStims{idx}.onsets    = trial_onsets(2);
    newStims{idx}.durations = trial_durations(2);
    newStims{idx}.amp       = 1;
    newStims{idx}.regressor_no_interest = 0;

    % --- Event 4: Cond 2 Train ---
    idx = 4;
    newStims{idx}.names     = [base_names{3} '_Train']; % e.g., 'HG_Train'
    newStims{idx}.onsets    = trial_onsets(3);
    newStims{idx}.durations = trial_durations(3);
    newStims{idx}.amp       = 1;
    newStims{idx}.regressor_no_interest = 0;

    % --- Event 5: Cond 2 Baseline ---
    idx = 5;
    newStims{idx}.names     = [base_names{3} '_Baseline']; % e.g., 'HG_Baseline'
    newStims{idx}.onsets    = baseline2_onset;
    newStims{idx}.durations = baseline2_duration;
    newStims{idx}.amp       = 1;
    newStims{idx}.regressor_no_interest = 0;

    % --- Event 6: Cond 2 Test ---
    idx = 6;
    newStims{idx}.names     = [base_names{4} '_Test']; % e.g., 'HG_Test'
    newStims{idx}.onsets    = trial_onsets(4);
    newStims{idx}.durations = trial_durations(4);
    newStims{idx}.amp       = 1;
    newStims{idx}.regressor_no_interest = 0;

    AllStims{i} = newStims;
end

% ---- Throw error if any exclusions occurred ----
exclude = unique(exclude);
exclude = exclude(exclude > 0);
if ~isempty(exclude)
    fprintf('--- DEBUG: exclusion details ---\n');
    for k = exclude
         trig_status = ~isempty(data_trigs{k}(1).Triggers) && numel(data_trigs{k}(1).Triggers) >= 8;
         time_status = ~isempty(raw(k).time) && isnumeric(raw(k).time);
         
         if ~trig_status
             msg = 'Reason: .tri file less than 8 triggers OR timestamp error.';
         elseif ~time_status
             msg = 'Reason: raw(i).time is empty or not numeric.';
         else
             msg = 'Reason: Unknown structural integrity issue.';
         end
         fprintf('Participant %d excluded. %s\n', k, msg);
    end
    error(['Participants missing data: ' num2str(exclude)]);
end

%% Working with nirs R
% r = nirs.core.Data;
% for i = 1:length(raw)
%     for j =1:4 % number of conditions
%         s = nirs.design.StimulusEvents;
%         s.name = AllStims{i}{j}.names;
%         s.onset = AllStims{i}{j}.onsets;
%         s.dur = AllStims{i}{j}.durations;
%         s.amp = AllStims{i}{j}.amp;
%         s.regressor_no_interest = AllStims{i}{j}.regressor_no_interest;
%         r(i).stimulus(s.name) = s;
%         % raw(i).stimulus(s.name) = s;
%         clear s
%     end
%     r(i).description    = raw(i).description;
%     r(i).data           = raw(i).data;
%     r(i).probe          = raw(i).probe;
%     r(i).time           = raw(i).time;
%     r(i).Fm             = raw(i).Fm;
%     r(i).auxillary      = raw(i).auxillary;
%     r(i).demographics   = raw(i).demographics;
% end
% Stim_Table = nirs.createStimulusTable(r); % create table of stim designs across subjects
% 
% %% Specify short-seperation channels
% for i = 1:length(r)
%     probe = r(i).probe;
%     probe.link.ShortSeperation = zeros(height(probe.link),1);
%     for j = 1:length(probe.link.detector)
%         if probe.link.detector(j) > 28
%             probe.link.ShortSeperation(j) = 1;
%         end
%         r(i).probe = probe;
%     end
% end
% 
% %% Conversions
% jobs         = nirs.modules.FixNaNs();
% jobs         = nirs.modules.OpticalDensity(jobs);
% jobs         = nirs.modules.Resample(jobs);
% jobs.Fs      = 1; % resample
% jobs         = nirs.modules.BeerLambertLaw(jobs);
% Hb = jobs.run(r);
% Hb(1).draw
% 
% %% Pre-processing & First-level analysis
% jobs         = nirs.modules.AddAuxRegressors();
% jobs.label   = ('aux');
% jobs         =nirs.modules.GLM(jobs);
% jobs.type    = 'AR-IRLS';
% jobs.AddShortSepRegressors = true;
% jobs         = nirs.modules.ExportData(jobs);
% jobs.Output  ='SubjStats';
% jobs.run(Hb)
% save('Hb_Final','SubjStats')
% 
% %% Second-level analysis
% j = nirs.modules.MixedEffects();
% j.formula = 'beta ~ -1 + cond + (1|Subject)'; %random effects
% j.dummyCoding = 'full';
% Group = j.run(SubjStats);
% 
% % (contrasts)
% Intrinsic = Group.ttest([0 0 1 -1]);
% Extrinsic = Group.ttest([0 1 0 -1]);
