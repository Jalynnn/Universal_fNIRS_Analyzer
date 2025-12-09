%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            SHINE LAB, CU BOULDER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the main analysis script for fNIRS neuroimaging data analysis.
% 2024 - Institute of Cognitive Science, CU-Boulder - James Crum  

% Add NIRS Toolbox to the Matlab path
clear, clc 
addpath(genpath('C:\Program Files\MATLAB\R2023a\toolbox\nirs-toolbox-master'));

% Initialize data directory
Directory = 'C:\Postdoc_CU\Navy\Data\raw'; 
cd(Directory(1:end-4));

% Load raw fNIRS intensities
rawDir = dir(Directory); rawDir = rawDir(3:end,:);
raw = nirs.io.loadDirectory(Directory, {'Subject'});
demographics = nirs.createDemographicsTable(raw);

% Load conditions
conditions = readtable('Conditions.csv');

%% Extract, organize, & insert stimulus designs
data_trigs = {};
trigDir = dir([Directory(1:end-4) filesep 'Trigs']); trigDir = trigDir(3:end,:);
for  i = 1:length(trigDir)
    cd([Directory(1:end-4) filesep 'Trigs' filesep trigDir(i).name])
    trigfile = dir(pwd); trigfile = trigfile(3:end,:);
    trigData = load(trigfile.name);
    data_trigs{i,1} = trigData(:,2);
end

% Remove particiants with insufficient timestamps across vars
a=1; idx = [];
for i = 1:length(data_trigs)
    [r c] = size(data_trigs{i});
    if r < 8
        idx(a) = i;
    end
    a=a+1;
end
idx(idx == 0) = [];

data_trigs(idx) = [];
raw(idx) = [];
conditions = table2cell(conditions);
conditions(idx,:) = [];

% Convert condition codes to string labels
for i = 1:length(conditions)
    for j = 2:5
        if str2num(conditions{i,j}(2)) == 2
            conditions{i,j} = 'LE_LI';
        end
        if str2num(conditions{i,j}(2)) == 6
            conditions{i,j} = 'HE_HI';
        end
        if str2num(conditions{i,j}(2)) == 3
            conditions{i,j} = 'LE_HI';
        end
        if str2num(conditions{i,j}(2)) == 4
            conditions{i,j} = 'HE_LI';
        end
    end
end

% Insert stimulus designs
a=1;
for i = 1:length(raw)
    starts = data_trigs{i}(1:2:end);
    ends   = data_trigs{i}(2:2:end);

    conds(1,1) = conditions(i,2);
    conds(2,1) = conditions(i,3);
    conds(3,1) = conditions(i,4);
    conds(4,1) = conditions(i,5);

    if length(raw(i).time) < starts(4)
        exclude(a) = i;
        a=a+1;
        continue
    end

    onsets(1,1) = raw(i).time(starts(1));
    onsets(2,1) = raw(i).time(starts(2));
    onsets(3,1) = raw(i).time(starts(3));
    onsets(4,1) = raw(i).time(starts(4));

    durations(1,1) = raw(i).time(ends(1)) - raw(i).time(starts(1));
    durations(2,1) = raw(i).time(ends(2)) - raw(i).time(starts(2));
    durations(3,1) = raw(i).time(ends(3)) - raw(i).time(starts(3));
    durations(4,1) = raw(i).time(ends(4)) - raw(i).time(starts(4));

    for j = 1:4 % conditions
        newStims{j}.names     = conds{j};
        newStims{j}.onsets    = onsets(j);
        newStims{j}.durations = durations(j);
        newStims{j}.amp       = ones(size(durations(j)));
        newStims{j}.regressor_no_interest = 0;
    end
    AllStims{i,1} = newStims;
end

% Exclude subjects with missing data
raw(exclude) = []; %nirs data
AllStims(exclude) = []; % behavioral data

r = nirs.core.Data;
for i = 1:length(raw)
    for j =1:4 % number of conditions
        s = nirs.design.StimulusEvents;
        s.name = AllStims{i}{j}.names;
        s.onset = AllStims{i}{j}.onsets;
        s.dur = AllStims{i}{j}.durations;
        s.amp = AllStims{i}{j}.amp;
        s.regressor_no_interest = AllStims{i}{j}.regressor_no_interest;
        r(i).stimulus(s.name) = s;
        % raw(i).stimulus(s.name) = s;
        clear s
    end
    r(i).description    = raw(i).description;
    r(i).data           = raw(i).data;
    r(i).probe          = raw(i).probe;
    r(i).time           = raw(i).time;
    r(i).Fm             = raw(i).Fm;
    r(i).auxillary      = raw(i).auxillary;
    r(i).demographics   = raw(i).demographics;
end
Stim_Table = nirs.createStimulusTable(r); % create table of stim designs across subjects

% Specify short-seperation channels
for i = 1:length(r)
    probe = r(i).probe;
    probe.link.ShortSeperation = zeros(height(probe.link),1);
    for j = 1:length(probe.link.detector)
        if probe.link.detector(j) > 28
            probe.link.ShortSeperation(j) = 1;
        end
        r(i).probe = probe;
    end
end

% Conversions
jobs         = nirs.modules.FixNaNs();
jobs         = nirs.modules.OpticalDensity(jobs);
jobs         = nirs.modules.Resample(jobs);
jobs.Fs      = 1; % resample
jobs         = nirs.modules.BeerLambertLaw(jobs);
Hb = jobs.run(r);
Hb(1).draw

% Pre-processing & First-level analysis
jobs         = nirs.modules.AddAuxRegressors();
jobs.label   = ('aux');
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

% (contrasts)
Intrinsic = Group.ttest([0 0 1 -1]);
Extrinsic = Group.ttest([0 1 0 -1]);

