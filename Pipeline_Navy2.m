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

    if cellSize < 10
        idx(a) = i;
        a = a + 1;
    end
end

% remove collected participants (if any)
if ~isempty(idx)
    error(['Error: The following participants have fewer than 10 triggers: ' ...
            num2str(idx)]);

    % % Which ptps to remove - Those with less than 10 triggers
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

    % Check for sufficient trigger points (must be >= 10: 4 trials * 2 + 2 baselines * 1)
    if isempty(trigVec) || numel(trigVec) < 10 || isempty(triTimesSec)
        exclude(a) = i; a = a+1; continue
    end

    % Secondary Exclusion Check (Raw NIRS Data Integrity)
    if isempty(raw(i).time) || ~isnumeric(raw(i).time)
        exclude(a) = i; a = a+1; continue
    end

    % 1. Separate Baseline Markers and Trial Triggers
    
    % The baseline *end* markers are triggers 31 and 32
    is_baseline_marker = (trigVec == 31) | (trigVec == 32); 
    
    % Baseline marker times (end of the 90s baseline period) and values
    baseline_marker_times = triTimesSec(is_baseline_marker);
    baseline_codes        = trigVec(is_baseline_marker);
    
    % Experimental trial trigger times and values
    is_trial = ~is_baseline_marker;
    trial_onsets_times = triTimesSec(is_trial);
    trial_codes        = trigVec(is_trial);

    % Sanity check: Should have 2 baseline markers and 8 trial triggers
    if length(baseline_marker_times) ~= 2 || length(trial_onsets_times) ~= 8
        warning('Participant %d: Expected 2 baseline and 8 trial triggers, found %d and %d.', ...
            i, length(baseline_marker_times), length(trial_onsets_times));
        exclude(a) = i; a = a+1; continue;
    end
    
    % 2. Extract Trial Onsets and Durations 
    % The experimental trials are still paired (odd = start, even = end)
    starts = 1:2:length(trial_codes);
    ends   = 2:2:length(trial_codes);

    % Onsets for the 4 experimental trials 
    trial_onsets    = trial_onsets_times(starts);
    trial_end_times = trial_onsets_times(ends); 
    trial_durations = trial_end_times - trial_onsets;

    % 3. Define Condition Labels
    % C1 is conditions{i,2} (e.g., 'LG'), C2 is conditions{i,3} (e.g., 'HG')
    cond_name_1 = conditions{i,2};
    cond_name_2 = conditions{i,3};

    % Use these base names to create the explicit stimulus names
    % Order: C1 Train, C1 Test, C2 Train, C2 Test
    base_names = { cond_name_1, cond_name_1, cond_name_2, cond_name_2 };

    % Need 6 total events: 4 trials + 2 baselines
    newStims = cell(6,1); 

    % 4. Calculate Baselines (90s BEFORE the marker) and Assign to Conditions

    % Determine which name corresponds to which code for the baseline
    if strcmp(cond_name_1, 'LG')
        lg_cond_name = cond_name_1;
        hg_cond_name = cond_name_2;
    else 
        lg_cond_name = cond_name_2;
        hg_cond_name = cond_name_1;
    end
    
    % Create the Baseline events structure
    baseline_events = struct('names', {}, 'onsets', {}, 'durations', {}, 'amp', {}, 'regressor_no_interest', {});
    
    for k = 1:2
        code = baseline_codes(k);
        marker_time = baseline_marker_times(k);
        
        % CRUCIAL CHANGE: Baseline starts 90 seconds BEFORE the marker time
        baseline_duration = 90;
        baseline_onset = marker_time - baseline_duration; 
        
        if code == 31 % LG Condition Baseline
            b_name = [lg_cond_name '_Baseline'];
        elseif code == 32 % HG Condition Baseline
            b_name = [hg_cond_name '_Baseline'];
        else
            error('Unexpected baseline code: %d', code);
        end
        
        baseline_events(k).names     = b_name; 
        baseline_events(k).onsets    = baseline_onset;
        baseline_events(k).durations = baseline_duration; % Fixed 90s duration
        baseline_events(k).amp       = 1;
        baseline_events(k).regressor_no_interest = 0;
    end
    
    % 5. Create Trial Events
    trial_events = cell(4, 1);
    for j = 1:4
        if j==1 || j==3; label='Train'; else; label='Test'; end
        trial_events{j}.names = [base_names{j} '_' label]; 
        
        trial_events{j}.onsets    = trial_onsets(j);
        trial_events{j}.durations = trial_durations(j);
        trial_events{j}.amp       = 1;
        trial_events{j}.regressor_no_interest = 0;
    end
    
    % 6. Combine and Sort All 6 Events Chronologically

    % Convert the 1x2 baseline struct array into a 2x1 cell array of structs
    baseline_events_cell = num2cell(baseline_events(:)); 
    % Note: `(:)` forces it into a column vector first, then `num2cell` wraps each struct in a cell.
    
    % Combine the two column cell arrays (2 baseline events + 4 trial events = 6 events)
    combined_events = [baseline_events_cell; trial_events]; % CORRECTED CONCATENATION

    % Extract all onsets from the combined events for sorting
    all_onsets_temp = cellfun(@(x) x.onsets, combined_events);

    % Sort the combined events based on their onset time
    [~, sorted_idx] = sort(all_onsets_temp);
    
    newStims = combined_events(sorted_idx);

    AllStims{i} = newStims;
end

% ---- Throw error if any exclusions occurred ----
exclude = unique(exclude);
exclude = exclude(exclude > 0);
if ~isempty(exclude)
    fprintf('--- DEBUG: exclusion details ---\n');
    for k = exclude
         trig_status = ~isempty(data_trigs{k}(1).Triggers) && numel(data_trigs{k}(1).Triggers) >= 10;
         time_status = ~isempty(raw(k).time) && isnumeric(raw(k).time);

         if ~trig_status
             msg = 'Reason: .tri file less than 10 triggers OR timestamp error.';
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
% Loop through each subject's raw data
for i = 1:length(raw)
    
    % --- CRITICAL FIX: CLEAR EXISTING STIMULUS DICTIONARY ---
    % The raw data loading function (nirs.io.loadDirectory) often inserts 
    % raw numeric triggers. We must clear these out before inserting the 
    % clean, named conditions from AllStims.
    % raw(i).stimulus = nirs.core.Dictionary; 
    % raw(i).stimulus = containers.Map;
    raw(i).stimulus = Dictionary;
    % --------------------------------------------------------
    
    % Access the cell array of stimulus structs for the current subject
    current_stims = AllStims{i};
    
    % Loop through all 6 stimulus events (2 Baselines + 4 Trials)
    for j = 1:length(current_stims) 
        
        % Create a new stimulus event object for the NIRS Toolbox
        s = nirs.design.StimulusEvents;
        
        % Check if the stimulus data is a struct and not empty
        stim_data = current_stims{j};
        if isstruct(stim_data)
            
            % Populate the stimulus event object 's' with the data you created
            s.onset = stim_data.onsets;
            s.dur = stim_data.durations;
            
            % --- FINAL NAME SANITIZATION CHECK ---
            temp_name = stim_data.names;
            
            % Sanitize any numeric name, just in case (though AllStims should be clean)
            if isstrprop(char(temp_name(1)), 'digit')
                 s.name = ['Condition_' temp_name];
                 warning('Sanitized numeric name in AllStims for subject %d, originally: %s', i, temp_name);
            else
                 s.name = temp_name; 
            end
            % ---------------------------------------
            
            % Force amplitude to a vector of ones (addressing the secondary concern)
            if ~isempty(s.onset)
                s.amp = ones(size(s.onset)); 
            else
                s.amp = 1;
            end
            
            % IMPORTANT: Insert the stimulus object into the raw data's stimulus dictionary
            raw(i).stimulus(s.name) = s; 
            
            clear s % Clear the temp variable before the next iteration
        else
            warning('Subject %d, Stimulus %d: Expected struct, but found unexpected type. Skipping.', i, j);
        end
    end
    
    % Optional: Add a check to ensure the stimulus field was populated
    if length(raw(i).stimulus.keys) ~= 6
        warning('Subject %d: Expected 6 stimuli, found %d after insertion.', i, length(raw(i).stimulus.keys));
    end
end

% Create a table of the stimulus designs across all subjects for validation
Stim_Table = nirs.createStimulusTable(raw); 
% Display the first few rows of the table to confirm the data looks correct
disp('*** Stimulus Table Head: ***');
disp(head(Stim_Table));

%% Specify short-seperation channels
for i = 1:length(raw) % <-- CORRECTED: Use 'raw' instead of the undefined 'r'
    probe = raw(i).probe;
    probe.link.ShortSeperation = zeros(height(probe.link),1);
    
    for j = 1:height(probe.link) % <-- CORRECTED: Use height(probe.link) for links
        if probe.link.detector(j) > 28
            probe.link.ShortSeperation(j) = 1;
        end
        raw(i).probe = probe;
    end
end
disp('*** SS Channels flagged on RAW data. ***');

%% Conversions
jobs         = nirs.modules.FixNaNs();
jobs         = nirs.modules.OpticalDensity(jobs);
jobs         = nirs.modules.Resample(jobs);
jobs.Fs      = 1; % resample
jobs         = nirs.modules.BeerLambertLaw(jobs);
Hb = jobs.run(raw);
Hb(1).draw

%% Pre-processing & First-level analysis
jobs         = nirs.modules.AddAuxRegressors();
% jobs.label   = ('aux');
jobs.label = {'aux'};
jobs         =nirs.modules.GLM(jobs);
jobs.type    = 'AR-IRLS';
jobs.AddShortSepRegressors = true;
jobs         = nirs.modules.ExportData(jobs);
jobs.Output  ='SubjStats';
jobs.run(Hb)
save('Hb_Final','SubjStats')

%% Second-level analysis
j = nirs.modules.MixedEffects();
j.formula = 'beta ~ -1 + cond + (1|Subject)'; %random effects
j.dummyCoding = 'full';
Group = j.run(SubjStats);

%% (contrasts)
% LG_Train; HG_Train; LG_Test; HG_Test
% White board

% HG_Train; HG_Baseline; HG_Test; LG_Train; LG_Baseline; LG_Test
c1 = Group.ttest([1 0 0 0 0 0]);
c2 = Group.ttest([0 0 1 0 0 1]);
c3 = Group.ttest([1 0 0 -1 0 0]);
c4 = Group.ttest([0 0 1 0 0 1]);
c5 = Group.ttest([1 0 -1 -1 0 -1]);