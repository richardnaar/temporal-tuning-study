%%  EEG at York %% Visit to York\Code\temporal-tuning-study
%% pwelch(), welch's method (Hann taper + overlaping windows) 
clc; close all; clear all; 

%% 
do = struct('other',                 0, 'plotting',             0,...
            'montage',               1, 'singleElectrodes',     0,...    
            'cleanData',             0, 'saveCsv',              0,...
            'saveCsv2',              0, 'saveCsvSingle',        1);

%% Paths

addpath('C:\Users\Richard Naar\Documents\Matlab toolboxes\ANTeepimport1.13\'); % 
dataDir = 'C:\Users\Richard Naar\Documents\dok\ssvep\Visit to York\EEG data\';  % directory for datafiles
dataDirClean = 'C:\Users\Richard Naar\Documents\dok\ssvep\Visit to York\EEG data\clean\';  % directory for datafiles
behDir = 'C:\Users\Richard Naar\Documents\dok\ssvep\Visit to York\EEG data\behavioural\';  % directory for datafiles

impdir =    [dataDir, 'raw\'];     
eeglabDir = ' C:\Program Files\MATLAB\R2014a\toolbox\eeglab13_6_5b\';                    % eeglab program directory
addpath(eeglabDir); addpath(impdir); addpath(behDir)

% % % Determine where your m-file's folder is.
% % folder = fileparts(which('pop_blinker.m')); 
% % % Add that folder plus all subfolders to the path.
% % addpath(genpath(folder));

implist = dir([impdir, '*.cnt']);          % make a fresh list of the contents of the raw files directory (input of this loop)
implistM = dir([dataDirClean, '*.mat']);   
%% other
KbName('UnifyKeyNames');
escapeKey = KbName('ESCAPE');   
KbCheck;
%% clean, find triggers, down-sample, find blinks and save

while do.cleanData
%% open the eeglab for a sec to get the paths as needed (only for pop_blinker)
eeglab; close
%%
for subi = 1:length(implist) % 
%% load data
EEG = pop_loadeep_v4([impdir, implist(subi).name]);
%EEG = eeg_checkset( EEG );

fprintf('loading participant: %s \n', implist(subi).name);
% eegplot(EEG.data, 'events', EEG.event)
%% channel locations
% EEG.chanlocs = ???
% there is also that option
%  EEG = pop_chanedit(EEG, 'load',{locfile 'filetype' 'autodetect'});     % Edit the channel locations structure of an EEGLAB dataset
%% find triggers
% import trigger information
[~,~,raw] = xlsread([dataDir, 'ViewPixx triggers - 3.xls']); % 3
fprintf('loading trigger info... \n');

%% find the triggers % 
trigNum = cell2mat(raw(2:5,2:4));
events = zeros(size(EEG.event));
labels = [];
trialCount = 0;
counter = 0;

for e = 1:size(EEG.event,2)
    if  size(EEG.event(e).type,2) > 4 
       events(1,e) = str2num(EEG.event(e).type(1)); 
       [R, C] = find(trigNum == str2num(EEG.event(e).type(1))); % 
    else 
       events(1,e) = str2num(EEG.event(e).type); 
       [R, C] = find(trigNum == str2num(EEG.event(e).type));

    end

    if strcmp((EEG.event(e).type), '185')
        label = 'ITI';
        EEG.event(e).type = label;
        trialCount = trialCount + 1;
    elseif strcmp((EEG.event(e).type), '81')
        label = 'RESCOR';
        EEG.event(e).type = label;
    elseif strcmp(EEG.event(e).type, '28')
        label = 'RESERR';
        EEG.event(e).type = label;
    elseif strcmp(EEG.event(e).type, '27')
        label = 'pause';
        EEG.event(e).type = label;
    elseif isempty(R)
        label = str2num(EEG.event(e).type);
        EEG.event(e).type = num2str(label);
    else
        row = char( raw(R+1,1) ); row = row(1:5);
        col = char( raw(1,C+1) ); col = col(1:4);
        label = [col,'_' row];
        label(regexp(label, '[ ]'))= []; 
        EEG.event(e).type = label; 
        counter = 0;
    end
    
    if e > 1 && size(EEG.event(e-1).type,2) > 2 % ~isempty(EEG.event(e-1).type) && 
        if ~isempty(strfind(EEG.event(e).type, EEG.event(e-1).type(1:3)))
        EEG.event(e).type = 'erroneusMarker';
        end
    end
    
 labels{1,e} = label;
end
fprintf('finding triggers... \n');

%% DOWN-SAMPLE 
if EEG.srate == 500
else
    EEG = pop_resample(EEG, 500);
end
%% find blinks with blinker

params = checkBlinkerDefaults(struct(), getBlinkerDefaults(EEG)); % extractin default parameters

% these parameters can be saved and loaded from the structure if needed,
% but I'll leave it like this for time being
params.stdThreshold = 1.5; % default
params.subjectID = implist(subi).name;
params.signalTypeIndicator = 'UseLabels'; % 'UseNumbers'; % 'UseLabels'
params.signalLabels = {'fp1'  'f1'  'fp2'  'fz'  'fpz'  'f3'  'f4'  'f2', 'veog'}; % find the best electrode from that set (added veog to the defaults)
params.blinkerSaveFile = [dataDir 'blinks\blinkDump'];
params.blinkerDumpDir = [dataDir 'blinks\blinkDump'];
params.dumpBlinkerStructures = true; % this just in case  
params.showMaxDistribution = true;
params.dumpBlinkImages = false; 
params.dumpBlinkPositions = false;
params.keepSignals = false;
params.keepSignals = false;

%save([dataDir, 'params-2AFC-TFUR-2019'], 'params', '-v7.3');

[~, com, blinks, blinkFits, blinkProperties, blinkStatistics, params] = pop_blinker(EEG, params); % struct()

% save the blink properties
save([[dataDir 'blinks\'] implist(subi).name '-bProp.mat'], 'blinkProperties', '-v7.3');
fprintf(['saveing blink properties as:' implist(subi).name '-bProp.mat \n']);

%eegplot(EEG.data, 'events', EEG.event)

%% (was moved to the next step)
eogs = {'VEOG', 'HEOG'}; 

eogi = find(ismember({EEG.chanlocs.labels}, eogs(2)));
EOG = EEG.data(eogi,:); EOGlocs = EEG.chanlocs(eogi);

EEG.data = EEG.data(1:end-2, :); EEG.chanlocs = EEG.chanlocs(1:end-2); EEG.nbchan = 64;

% EEG = clean_drifts(EEG, [0.25 0.75]); % Removes drifts from the data using a forward-backward high-pass filter.

% Maximum tolerated flatline duration. In seconds. If a channel has a longer flatline than this, it will be considered abnormal. 
EEG = clean_flatlines(EEG, 5);

[EEG, rej] = pop_rejchan(EEG, 'elec',[1:EEG.nbchan] ,'threshold',8,'norm','on','measure','kurt');

% rejElecs{subi} = rej;

%% compute average reference
% eegplot(EEG.data, 'events', EEG.event) 
EEG_ref = mean(EEG.data(1:end,:),1); % 1:64

%% adding EOG channel back
EEG.data = [EEG.data; EOG]; EEG.chanlocs(end+1) = EOGlocs; EEG.nbchan = EEG.nbchan + 1;

EEG.data = bsxfun(@minus, EEG.data(1:end,:), squeeze(EEG_ref)); % average reference (excluding EOG from the mean)

fprintf(['Averge reference (excluding EOG from the mean) \n'])

%% epoch with eeglab
% eventType = {{'STIM_HIGH' 'STIM_LOW' 'STIM_RND_H' 'STIM_RND_L'} [-7 0] [-200 0]}; 
% eventType = {{'STIM_HIGH'} [-7 0] [-200 0]}; 

% EEG = pop_epoch(EEG, eventType{1}, eventType{2}, 'epochinfo', 'yes');  % extract epochs 
% 
% EEG = pop_rmbase(EEG, []);                                    % (useing whole epoch as baseline) does it make sense to remove baseline eventType{3}
%% convert to FieldTrip format
% kuku = eeglab2fieldtrip(EEG,'preprocessing');
%% save the data

save([dataDirClean, implist(subi).name(1:end-4), '-clean.mat'], 'EEG', '-v7.3');
fprintf(['Saveing the dataset as: ' implist(subi).name(1:end-4), '-clean.mat \n'])


end
do.cleanData = 0;
end

while do.montage
% eegplot(EEG.data, 'events', EEG.event)
%% epoch and save single trial and grand average SNR and amplitude data, large data structure (p) and ress

doRess = 1;
do_p_struct = 0;
doSingleTrial = 1;
plotPlinks = 0;
findBlinks = 0;
stimLocked = 0; % stim locked events allow to preserve data near the stim event

if do_p_struct; doRess = 1; doSingleTrial = 1; doAvgAmp = 1; doAvgSNR = 1; end
%% some electrode sets

%left_frontal = {'Fp1', 'AF3', 'AF7', 'F1', 'F3', 'F5', 'F7'};
% right_frontal = {'Fp2', 'AF4', 'AF8', 'F2', 'F4', 'F6', 'F8'}; 
% left_fronto_central = {'FC1', 'FC3', 'FC5', 'C1', 'C3', 'C5', 'CP1', 'CP3', 'CP5'};
% right_fronto_central = {'FC2', 'FC4', 'FC6', 'C2', 'C4', 'C6', 'CP2', 'CP4', 'CP6'}; 
% left_temporal = {'FT7', 'T7', 'TP7', 'TP9'};
%right_temporal = {'FT8', 'T8', 'TP8', 'TP10'};
% left_parieto_occipital = {'P1', 'P3', 'P5', 'P7', 'PO3', 'PO7', 'O1'};
% right_parieto_occipital = {'P2', 'P4', 'P6', 'P8', 'PO4', 'PO8', 'O2'};
% fronto_central = {'Fpz', 'Fz'};
% central = {'FCz', 'Cz', 'CPz'}; 
occipito_central = {'Pz', 'POz', 'Oz', 'O1', 'O2'}; 
% occipito_central2 = {'Pz', 'POz', 'Oz', 'O1', 'O2'}; 
% right_fronto_temporal = {'FT8', 'T8', 'TP8', 'TP10'}; 
mt_left = {'P5', 'P7', 'PO7'}; % BA 19, 37 39 http://www.brainm.com/software/pubs/dg/BA_10-20_ROI_Talairach/nearesteeg.htm
% mt_right = {'P6', 'P8', 'PO4'};
% mt_right = {'P6', 'P8', 'PO4'};
mt_right2 = {'P6', 'P8', 'PO8'};
V3 = {'PO4', 'PO8'};

current_set = V3;

%% loop through all the datasets in implistM (pre-processed data) one by one
manualCheck = 0; % check manually for artefacts
for subi = 1:length(implistM) 
% empty them after each iteration
meanComplexFFT = [];
thisEventBlinks = {};
    
%% load psychophysics data (this takes currently an assumption that the EEG failes are in the same order as csv files)

files = dir([behDir , '*.csv']);
filename = fullfile(behDir,files(subi).name); % find the corresponding file

behDat = readtable(filename,...
    'Delimiter',',','ReadVariableNames',false); % load the data

names = table2array(behDat(1,4:41)); % take the variable names from the first row
names = strrep(names, '.',''); names = strrep(names, '_',''); names = strrep(names, ' ','');

behDat = behDat(3:end,4:41);
behDat.Properties.VariableNames = names; % rename the variables
behDatConds = {'high', 'low', 'high50', 'low50'};
%% load EEG
EEG = []; 

load([dataDirClean, implistM(subi).name]); % load data

fprintf('loading participant: %s \n', implistM(subi).name);

srate = EEG.srate;  % sampling frequency 
EEG_raw = EEG;
%% filter 
lowcutoff = 1;
revfilt = 0; % [0|1] reverse filter (i.e. bandpass filter to notch filter). {default 0}
plotfreqz = 0;
minphase = false; % defalut

[EEG, com, b] = pop_eegfiltnew(EEG, 1, 125, [], revfilt, [], plotfreqz, minphase);
% eegplot(EEG.data, 'events', EEG.event) % see the raw file

%% extracting the blink structure
if findBlinks
    blinkList = dir([dataDir 'blinks\', '*.mat']); 
    % load([dataDir 'blinks\' blinkList(end).name])
    % fprintf('loading dataset: %s', blinkList(end).name); fprintf('\n') 

    bFileCounter = 0; fileNotFound = 1;
    while fileNotFound
    bFileCounter = bFileCounter + 1;
        if strcmp( blinkList( bFileCounter ).name(1:7) , implistM(subi).name(1:7) ) ~= 1 && bFileCounter <= length(blinkList)
            % do nothing
        elseif strcmp( blinkList(bFileCounter).name(1:7) , implistM(subi).name(1:7) ) == 1
            load([dataDir 'blinks\' blinkList(bFileCounter).name])
            fprintf('loading dataset: %s', blinkList(bFileCounter).name); fprintf('\n') 
            fileNotFound = 0;
        else
            fprintf('Blinker file not found')
        end
    end
end
%% epoch parameters


if stimLocked == 1
    event = {'STIM_HIGH' 'STIM_LOW' 'STIM_RND_H' 'STIM_RND_L'}; 
else
    event = {'PRED_HIGH' 'PRED_LOW' 'PRED_RND_H' 'PRED_RND_L'};
end

% other events
% event = {{'STIM_HIGH', 'STIM_LOW'}, {'STIM_RND_L', 'STIM_RND_H'}}; 
% event = {'PRED_HIGH' 'PRED_LOW' 'PRED_RND_H' 'PRED_RND_L'}; 
% event = {'CUE_HIGH' 'CUE_LOW' 'CUE_RND_H' 'CUE_RND_L'};
% event = {'ITI'};

% this is the transient response that will be subtracted from the epoched data
% note that in majority of cases it leads to noninteger durations
transient = 0.5; % this is in seconds

% subtract transient from 10 delay periods used in the study (in seconds)
% and multiply it by the sampling rate
if stimLocked == 1
    eventDurs = ([1.5, 2, 2.5, 3]-transient ) * srate;
else
    eventDurs = ( [ 2, 2.5, 3, 3.5, 4, 5, 6, 7, 8 ]-transient ) * srate;
end

% find maximum duration in samples that will be divisible by sampling rate
% (this is needed for coherence averaging)
eventDurs = eventDurs - mod(eventDurs, srate);

% number of trials to skip in each event category
fromTrial = 0; 

%% save same data for the p structure
if do_p_struct
    p.subject(subi).conditionLabels  = event;
    % p.conditionLabels = event;
    p.transient = transient;
    p.fromTrial = fromTrial;
    p.subject(subi).psychoLabels = behDatConds; 
    p.subject(subi).Id = implistM(subi).name;
end
%%
for eventIndx = 1:length(event) % loop through the event categories

% Find all the events corresponding to the eventIndex
thisEvent = event{eventIndx}; % pick events from predefined event categories
nEvents=length(EEG.event); % get the number of all the events in the dataset
currentIndex=1; % set the eventindex to be one before entering the loop 
% (this will be increased by one each time we find the predefined event from 
% the eventlist, after exclusion of the number of trials set in freomTrial)
eventCounter = 0; % set the counter to exclude first fwe trials

for thisEventIndx=1:nEvents % loop through all the events
    if  strfind(EEG.event(thisEventIndx).type, char(thisEvent)) % any(ismember(thisEvent, EEG.event(thisEventIndx).type)) % to compare cued to non-cued
        eventCounter = eventCounter + 1;
        if eventCounter > fromTrial %
            if stimLocked == 1
                eventOffset(currentIndex) = EEG.event(thisEventIndx).latency; % find event offset
                eventOnset(currentIndex) = EEG.event(thisEventIndx-1).latency; % find previous event oncet (start of the epoch)
            else
                eventOnset(currentIndex) = EEG.event(thisEventIndx).latency; % find event oncet
                eventOffset(currentIndex) = EEG.event(thisEventIndx+1).latency; % find next event oncet (end of the epoch)
            end
            currentIndex=currentIndex+1;
        end
    end
end
 
fprintf('\nFound %d events \n',currentIndex-1);
%% select channels
% electrodes = {'O1', 'Oz', 'O2', 'HEOG'}; % HEOG will be removed latere (NB! hard coded) 'POz', 'P1', 'P2', 'PO5', 'PO3', 'PO4', 'PO6',

electrodes = {current_set{:}, 'HEOG'}; % HEOG will be removed latere (NB! hard coded) 
if do_p_struct; p.electrodes = electrodes(1:end-1); end % save to p

elec2plot = find(ismember({EEG.chanlocs.labels}, electrodes)); % find electrode indexes
elec2plotNames = {EEG.chanlocs(elec2plot).labels}; % save the electrode names to the variable

fprintf('\nNumber of electrodes aggregated (NB! last channel, which is EOG, will be discarded - hard coded):  %d ', length(elec2plot)); fprintf('\n')

% elec2plotNames 
%% prepare psychophysics data
condi = find(ismember(behDat.trials2label, behDatConds(eventIndx))); % find trial indexes for that event type
% behDat.trials2response(condi())
%% epoch the data
% ressProov = double.empty(4, 64, 3500,0);

% these will be over writen for each event type
ressData = {}; ressDataAvg = {};
thisCueEvent = 0;
thisEventBlinks = {};
meanComplexFFT = [];
tIndex = zeros(65,1);% []; % valid coherence trials
validRess = zeros(65,1); % valid ress trials
rtPsychoTrials = []; rtPsychoTrials = []; correctPsychoTrials = []; % not rally necessary but ok 
pauseForInspection = 0;
%%
while (thisCueEvent <= currentIndex-2) % for thisCueEvent = 1:(currentIndex-1)
     thisCueEvent = thisCueEvent + 1;
     
     startSample = round(eventOnset(thisCueEvent)+transient*srate); % *EEG.srate/1000 % Compensate for the fact that latencies are in ms but data are in samples. Is that really the case?   
     % +(EEG.srate/2)
     endSample = round(eventOffset(thisCueEvent)); % 
     
     trialDur = eventDurs(dsearchn(eventDurs', endSample-startSample)); % find the trial duration matching the current trial

     % redifine the end or the start of the epoch to be suitable for
     % coherent averaging 

     startSample = endSample-trialDur; % find the event oncet if the epoch is calculated from stim event       

          
    % Pull out a timeseries that follows/preceeds the event 

    allSamples = EEG.data(elec2plot,startSample:endSample-1); % channels
    
%% prepare ress data matrix
    if doRess
        allSamplesRess = EEG_raw.data(1:end-1,startSample:endSample-1); % channels % not including the ocular channel (end-1) 
        allSamplesRess = bsxfun(@minus, allSamplesRess, mean(EEG_raw.data(1:end-1,startSample-srate/2:startSample),2)); % 
        rebinnedDataRess = reshape(allSamplesRess, size(allSamplesRess,1), EEG.srate,trialDur/EEG.srate);    
        artefacts = ones(size(rebinnedDataRess,3),1); % accept all for RESS
    end
%     ressThisEvent = rebinnedDataRe; % not including the ocular channel (end-1) 

%     [artefacts, ~] = findManual(rebinnedDataRess, manualCheck, size(rebinnedDataRess,3), escapeKey, 1);


    if doRess && sum(artefacts) ~= 0
        validRess(thisCueEvent) = 1;
        ressData{thisCueEvent} = rebinnedDataRess(:,:,find(artefacts));
        ressDataAvg{thisCueEvent} = mean(rebinnedDataRess(:,:,find(artefacts)),3);
    end 
    
%% just to inspect the single trials    

rt = str2num(behDat.trueRT{condi(thisCueEvent+fromTrial)});
isCorrect = behDat.trials2response(condi(thisCueEvent+fromTrial));
contrast = str2num(behDat.trials2intensity{condi(thisCueEvent+fromTrial)});

subtrMean = bsxfun(@minus, allSamples, mean(allSamples,2)); % 
zscored = bsxfun(@rdivide, subtrMean, std(allSamples,0,2)); % 

%%
if plotPlinks && findBlinks
    figure(1)
    subplot(2,1,1)
    plot(allSamples')
    hold on
end

%% find and draw blinks for this trial
if findBlinks 
startSampleInSeconds = startSample/srate;
endSampleInSeconds = endSample/srate;


% extend the window size
extS = 0.01;     
extE = 0.01;
if (startSampleInSeconds  - extS) > 0
    startSampleInSeconds = startSampleInSeconds - extS;
end
    endSampleInSeconds = endSampleInSeconds + extE;
    
    
numElements = arrayfun(@(x) (x.peakTimeBlink >= startSampleInSeconds && x.peakTimeBlink <= endSampleInSeconds) , blinkProperties);
bIndx = find(numElements);

bcount = size(bIndx,2);

nrSegments = trialDur/srate;
notblinks = ones(nrSegments,1);
blinkPeaks = zeros(size(bIndx,2),1);
for ii = 1:size(bIndx,2)
    peakT = round( blinkProperties(bIndx(ii)).peakTimeBlink*500 ); % find current blink indexes
    dur = round( (blinkProperties(bIndx(ii)).durationBase*500) /2 ); % find the blinks left and right shoulder
    
    % recalibrate for the plot
     peakT = peakT - startSample; 
     blinkPeaks(ii,1) = peakT;
if plotPlinks    
    yl = ylim; 
%         for hpi = 1:size(blinkArtifacts,2) 
        hp = patch([peakT-dur peakT-dur peakT+dur peakT+dur ]...
            ,[yl(1) yl(2) yl(2) yl(1)],'k',...
        'facecolor', [.75,.75,.75],'edgecolor',[.75,.075,.275],...
        'erasemode','xor') ;
%         end
%     end 




    title(['subject: ',num2str(subi),'; trial: ', num2str(thisCueEvent+fromTrial),...
        '; condition: ', strrep(event{eventIndx}, '_',' '), '; is correct: ',...
        isCorrect{:}, '; RT: ', num2str(rt), ' blinks: ', num2str(bcount)]) % 
%     set(gca,'xlim', [0 3500], 'FontSize',12)
end

if pauseForInspection && plotPlinks
    pause()
end
 

%% find a vector containing boolean values for segment rejection (1 keep the trial)
if ii == size(bIndx,2)
    
a = ceil( blinkPeaks/srate);
a(a == 0) = []; a(a > nrSegments) = [];

seg = zeros(nrSegments,1);     
seg(a,1) = 1;

notblinks = seg < 1;

end
% txt = find(notblinks == 0)

end

end
if plotPlinks && findBlinks; hold off; end

%% extract some data
if isempty(rt)
    rtPsychoTrials(thisCueEvent+fromTrial) = nan;
else
    rtPsychoTrials(thisCueEvent+fromTrial) = rt;
end

contrastPsychoTrials(thisCueEvent+fromTrial) = contrast;

correctPsychoTrials(thisCueEvent+fromTrial) = str2num(isCorrect{:});

%% prepare fft data   
    if size(elec2plot,2)-1 > 1
        allSamples=squeeze(mean(allSamples(1:end-1, :))); % mean of channels (not including ocular channel)
    else
        allSamples=squeeze(allSamples(1:end-1, :)); 
    end
    
    % We can now bin into 1 second intervals
    rebinnedData=reshape(allSamples, EEG.srate,trialDur/EEG.srate);
    fftRebinned=fft(rebinnedData); % Perform FFT down time
    % Now the magic happens :) We can average coherently across bins to remove
    % noise
%%
goBack = 0;
%[~, goBack] = findManual(rebinnedData, manualCheck, trialDur/EEG.srate, escapeKey, 2);
% notblinks is a vector where 0 indicates a blink

if findBlinks == 0
    notblinks = 1:size(rebinnedData,2);
end

if sum(notblinks) ~= 0 % If there is at least one segment without a blink
    rebinnedData = rebinnedData(:,find(notblinks));
    
%     blinkIndex = blinkIndex + 1;
%     trialId(thisCueEvent) = thisCueEvent+fromTrial;
%     meanComplexFFTallTrials(:,blinkIndex) = mean(fftRebinned,2); % keep this for all the trials

    meanComplexFFT(:,thisCueEvent) = mean(fftRebinned,2); % nb abs % This is the mean complex FFT for this trial (averaged across bins)
    tIndex(thisCueEvent+fromTrial) = 1;


else  
    fprintf('no good segments found in this trial \n')
end

thisEventBlinks{thisCueEvent,1} = notblinks;
% 

if (thisCueEvent + goBack) > 0
thisCueEvent = thisCueEvent + goBack;
end

%meanComplexFFT(:,thisCueEvent) = mean(fftRebinned,2); % nb abs % This is the mean complex FFT for this trial (averaged across bins)
%thisEventBlinks{thisCueEvent,1} = blinks;   
    
end

% grandAverage{eventIndx}=mean((meanComplexFFT),2); % Average across trials - the abs means we do NOT still keep coherent information
grandAverage{subi, eventIndx} = mean(meanComplexFFT(:,find(tIndex)),2); % Average across trials - the abs means we do NOT still keep coherent information
if ~isempty(meanComplexFFT) && doSingleTrial
    allTheSingleTrials{subi, eventIndx} =  abs(meanComplexFFT(1:60,find(tIndex)) ).^2; % Here we lose coherence.
end

if findBlinks
    allBlinks{subi, eventIndx} = thisEventBlinks;
end

if do_p_struct
    % for the p
    p.subject(subi).condition(eventIndx).data.EEG.artefacts = thisEventBlinks; % allBlinks; 
    p.subject(subi).condition(eventIndx).data.EEG.validTrials = {find(tIndex)}; 
    p.subject(subi).condition(eventIndx).data.EEG.single.amplitude = allTheSingleTrials; 
    p.subject(subi).condition(eventIndx).data.EEG.average.amplitude = grandAverage; 

    p.subject(subi).condition(eventIndx).data.EEG.ressDataSegments = ressData;% squeeze( ressData(subi, eventIndx, find(validRess)) ); 
    p.subject(subi).condition(eventIndx).data.EEG.ressDataTrials = ressDataAvg;
    p.ressDataStructure = {'participant', 'condition', 'trials', 'time'};
    p.allTheSingleTrialsStructure = {'participant', 'condition', 'frequency', 'trials'};
    
    % rt, intesity, correct
    p.subject(subi).condition(eventIndx).data.psychophysics.correct = {correctPsychoTrials}; 
    p.subject(subi).condition(eventIndx).data.psychophysics.rt = {rtPsychoTrials};
    p.subject(subi).condition(eventIndx).data.psychophysics.contrast = {contrastPsychoTrials};
end

if doRess
    ressAllSegments.subject(subi).condition(eventIndx).data = ressData;% squeeze( ressData(subi, eventIndx, find(validRess)) ); 
    ressAllTrials.subject(subi).condition(eventIndx).data = ressDataAvg;% 
end

%% SNR (need to check this: 27.01.22)

hz = abs(grandAverage{subi, eventIndx}(2:44)).^2;
snrE = zeros(1,size(hz,1));
skipbins =  1; % 1 Hz, hard-coded! (not skipping)
numbins  = 2; %  2 Hz, also hard-coded!

% loop over frequencies and compute SNR
for hzi=( numbins+1 ):( length(hz)-numbins-1 )
    numer = hz(hzi);
    denom = rms( hz([( hzi-numbins ):(hzi-skipbins) (hzi+skipbins):(hzi+numbins)]) ); 
    snrE(hzi) =  (numer-denom)/(numer+denom); %numer./denom;
end  
%%
% %     subplot(2,2,eventIndx); hold on
% %     title(strrep(event{eventIndx}, '_',' '))

%    bar(abs(grandAverage(2:60)));

%     bar(abs(grandAverage{eventIndx}(2:60)));
allSnrE{subi,eventIndx} = snrE;

if do_p_struct; p.subject(subi).condition(eventIndx).data.EEG.average.snr = snrE; end

% %     bar(allSnrE{subi,eventIndx});

%     stem(abs(grandAverage{eventIndx}(2:60)),'k','linew',3,'markersize',2.5,'markerfacecolor','r')
%     bar(squeeze( mean(snrE(:, eventIndx, 1:end),1)));
end % next event
if do_p_struct; p.subject(subi).psychoPyTable = behDat; end
% pid = p.subject(subi);
% save([[dataDir 'individual\'] date p.subject(subi).Id '.mat'], 'pid', '-v7.3');

end % next subject

% save data ressData
if stimLocked == 1
    save([dataDir date '-TFUR-SNRs-V3-stimLocked.mat'], 'allSnrE', '-v7.3');
    if doSingleTrial; save([dataDir date '-TFUR-amp-V3-stimLocked.mat'], 'allTheSingleTrials', '-v7.3'); end
    save([dataDir date '-TFUR-grandAverage-V3-stimLocked.mat'], 'grandAverage', '-v7.3');
else
    save([dataDir date '-TFUR-SNRs-cueLocked-V3.mat'], 'allSnrE', '-v7.3');
    save([dataDir date '-TFUR-grandAverage-cueLocked-V3.mat'], 'grandAverage', '-v7.3');
    if doSingleTrial; save([dataDir date '-TFUR-amp-cueLocked-V3.mat'], 'allTheSingleTrials', '-v7.3'); end
end

if doRess
    save([dataDir date '-TFUR-ressAllSegments-cueLocked-V3.mat'], 'ressAllSegments', '-v7.3');
    save([dataDir date '-TFUR-ressAllTrials-cueLocked-V3.mat'], 'ressAllTrials', '-v7.3');
end

if do_p_struct
    save([dataDir date '-allData.mat'], 'p', '-v7.3');  
end
do.montage = 0;
end

while do.singleElectrodes
%% all electrode names    
all64 = {'Fp1','Fpz','Fp2','F7','F3','Fz','F4','F8','FC5','FC1','FC2','FC6','M1','T7','C3','Cz','C4','T8','M2','CP5','CP1','CP2'...
    'CP6','P7','P3','Pz','P4','P8','POz','O1','Oz','O2','AF7','AF3','AF4','AF8','F5','F1','F2','F6','FC3','FCz','FC4','C5'...
'C1','C2','C6','CP3','CPz','CP4','P5','P1','P2','P6','PO5','PO3','PO4','PO6','FT7','FT8','TP7','TP8','PO7','PO8'};
elecStruct = struct();
%% loop through all the datasets in implistM (pre-processed data) one by one
manualCheck = 0; % check manually for artefacts
for subi = 1:length(implistM) 
% empty them after each iteration
plotPlinks = 0;
meanComplexFFT = [];
thisEventBlinks = {};
    
%% load psychophysics data (this takes currently an assumption that the EEG failes are in the same order as csv files)

files = dir([behDir , '*.csv']);
filename = fullfile(behDir,files(subi).name); % find the corresponding file

behDat = readtable(filename,...
    'Delimiter',',','ReadVariableNames',false); % load the data

names = table2array(behDat(1,4:41)); % take the variable names from the first row
names = strrep(names, '.',''); names = strrep(names, '_',''); names = strrep(names, ' ','');

behDat = behDat(3:end,4:41);
behDat.Properties.VariableNames = names; % rename the variables
behDatConds = {'high', 'low', 'high50', 'low50'};
%% load EEG
EEG = []; 

load([dataDirClean, implistM(subi).name]); % load data

fprintf('loading participant: %s \n', implistM(subi).name);

srate = EEG.srate;  % sampling frequency 
EEG_raw = EEG;
%% filter 
lowcutoff = 1;
revfilt = 0; % [0|1] reverse filter (i.e. bandpass filter to notch filter). {default 0}
plotfreqz = 0;
minphase = false; % defalut

[EEG, com, b] = pop_eegfiltnew(EEG, 1, 125, [], revfilt, [], plotfreqz, minphase);
% eegplot(EEG.data, 'events', EEG.event) % see the raw file

%% extracting the blink structure
% 
blinkList = dir([dataDir 'blinks\', '*.mat']); 
% load([dataDir 'blinks\' blinkList(end).name])
% fprintf('loading dataset: %s', blinkList(end).name); fprintf('\n') 

bFileCounter = 0; fileNotFound = 1;
while fileNotFound
bFileCounter = bFileCounter + 1;
    if strcmp( blinkList( bFileCounter ).name(1:7) , implistM(subi).name(1:7) ) ~= 1 && bFileCounter <= length(blinkList)
        % do nothing
    elseif strcmp( blinkList(bFileCounter).name(1:7) , implistM(subi).name(1:7) ) == 1
        load([dataDir 'blinks\' blinkList(bFileCounter).name])
        fprintf('loading dataset: %s', blinkList(bFileCounter).name); fprintf('\n') 
        fileNotFound = 0;
    else
        fprintf('Blinker file not found')
    end
end

%% epoch parameters

stimLocked = 0; % stim locked events allow to preserve data near the stim event

if stimLocked == 1
    event = {'STIM_HIGH' 'STIM_LOW' 'STIM_RND_H' 'STIM_RND_L'}; 
else
    event = {'PRED_HIGH' 'PRED_LOW' 'PRED_RND_H' 'PRED_RND_L'};
end


% this is the transient response that will be subtracted from the epoched data
% note that in majority of cases it leads to noninteger durations
transient = 0.5; % this is in seconds

% subtract transient from 10 delay periods used in the study (in seconds)
% and multiply it by the sampling rate
if stimLocked == 1
    eventDurs = ([1.5, 2, 2.5, 3]-transient ) * srate;
else
    eventDurs = ( [ 2, 2.5, 3, 3.5, 4, 5, 6, 7, 8 ]-transient ) * srate;
end

% find maximum duration in samples that will be divisible by sampling rate
% (this is needed for coherence averaging)
eventDurs = eventDurs - mod(eventDurs, srate);

% number of trials to skip in each event category
fromTrial = 0; 

%% save same data for the p structure

%%
for eleci = 1:64
current_elctrode = all64(eleci);    
for eventIndx = 1:length(event) % loop through the event categories

% Find all the events corresponding to the eventIndex
thisEvent = event{eventIndx}; % pick events from predefined event categories
nEvents=length(EEG.event); % get the number of all the events in the dataset
currentIndex=1; % set the eventindex to be one before entering the loop 
% (this will be increased by one each time we find the predefined event from 
% the eventlist, after exclusion of the number of trials set in freomTrial)
eventCounter = 0; % set the counter to exclude first fwe trials

for thisEventIndx=1:nEvents % loop through all the events
    if  strfind(EEG.event(thisEventIndx).type, char(thisEvent)) % any(ismember(thisEvent, EEG.event(thisEventIndx).type)) % to compare cued to non-cued
        eventCounter = eventCounter + 1;
        if eventCounter > fromTrial %
            if stimLocked == 1
                eventOffset(currentIndex) = EEG.event(thisEventIndx).latency; % find event offset
                eventOnset(currentIndex) = EEG.event(thisEventIndx-1).latency; % find previous event oncet (start of the epoch)
            else
                eventOnset(currentIndex) = EEG.event(thisEventIndx).latency; % find event oncet
                eventOffset(currentIndex) = EEG.event(thisEventIndx+1).latency; % find next event oncet (end of the epoch)
            end
            currentIndex=currentIndex+1;
        end
    end
end
 
fprintf('\nFound %d events \n',currentIndex-1);
%% select channels
% electrodes = {'O1', 'Oz', 'O2', 'HEOG'}; % HEOG will be removed latere (NB! hard coded) 'POz', 'P1', 'P2', 'PO5', 'PO3', 'PO4', 'PO6',

electrodes = {current_elctrode{:}, 'HEOG'};

elec2plot = find(ismember({EEG.chanlocs.labels}, electrodes)); % find electrode indexes
elec2plotNames = {EEG.chanlocs(elec2plot).labels}; % save the electrode names to the variable

fprintf('\nNumber of electrodes aggregated (NB! last channel, which is EOG, will be discarded - hard coded):  %d ', length(elec2plot)); fprintf('\n')

if length(elec2plot) > 1
% elec2plotNames 
%% prepare psychophysics data
condi = find(ismember(behDat.trials2label, behDatConds(eventIndx))); % find trial indexes for that event type
% behDat.trials2response(condi())
%% epoch the data
% ressProov = double.empty(4, 64, 3500,0);

% these will be over writen for each event type
ressData = {}; ressDataAvg = {};
thisCueEvent = 0;
thisEventBlinks = {};
meanComplexFFT = [];
tIndex = zeros(65,1);% []; % valid coherence trials
validRess = zeros(65,1); % valid ress trials
rtPsychoTrials = []; rtPsychoTrials = []; correctPsychoTrials = []; % not rally necessary but ok 
pauseForInspection = 0;
%%
while (thisCueEvent <= currentIndex-2) % for thisCueEvent = 1:(currentIndex-1)
     thisCueEvent = thisCueEvent + 1;
     
     startSample = round(eventOnset(thisCueEvent)+transient*srate); % *EEG.srate/1000 % Compensate for the fact that latencies are in ms but data are in samples. Is that really the case?   
     % +(EEG.srate/2)
     endSample = round(eventOffset(thisCueEvent)); % 
     
     trialDur = eventDurs(dsearchn(eventDurs', endSample-startSample)); % find the trial duration matching the current trial

     % redifine the end or the start of the epoch to be suitable for
     % coherent averaging 

     startSample = endSample-trialDur; % find the event oncet if the epoch is calculated from stim event       

          
    % Pull out a timeseries that follows/preceeds the event 

    allSamples = EEG.data(elec2plot,startSample:endSample-1); % channels
    
%% just to inspect the single trials    

rt = str2num(behDat.trueRT{condi(thisCueEvent+fromTrial)});
isCorrect = behDat.trials2response(condi(thisCueEvent+fromTrial));
contrast = str2num(behDat.trials2intensity{condi(thisCueEvent+fromTrial)});

subtrMean = bsxfun(@minus, allSamples, mean(allSamples,2)); % 
zscored = bsxfun(@rdivide, subtrMean, std(allSamples,0,2)); % 

%%
if plotPlinks
    figure(1)
    subplot(2,1,1)
    plot(allSamples')
    hold on
end
%% find and draw blinks for this trial
 
startSampleInSeconds = startSample/srate;
endSampleInSeconds = endSample/srate;


% extend the window size
extS = 0.01;     
extE = 0.01;
if (startSampleInSeconds  - extS) > 0
    startSampleInSeconds = startSampleInSeconds - extS;
end
    endSampleInSeconds = endSampleInSeconds + extE;
    
    
numElements = arrayfun(@(x) (x.peakTimeBlink >= startSampleInSeconds && x.peakTimeBlink <= endSampleInSeconds) , blinkProperties);
bIndx = find(numElements);

bcount = size(bIndx,2);

nrSegments = trialDur/srate;
notblinks = ones(nrSegments,1);
blinkPeaks = zeros(size(bIndx,2),1);
for ii = 1:size(bIndx,2)
    peakT = round( blinkProperties(bIndx(ii)).peakTimeBlink*500 ); % find current blink indexes
    dur = round( (blinkProperties(bIndx(ii)).durationBase*500) /2 ); % find the blinks left and right shoulder
    
    % recalibrate for the plot
     peakT = peakT - startSample; 
     blinkPeaks(ii,1) = peakT;
     
if plotPlinks    
    yl = ylim; 

    hp = patch([peakT-dur peakT-dur peakT+dur peakT+dur ]...
    ,[yl(1) yl(2) yl(2) yl(1)],'k',...
    'facecolor', [.75,.75,.75],'edgecolor',[.75,.075,.275],...
    'erasemode','xor') ;

    title(['subject: ',num2str(subi),'; trial: ', num2str(thisCueEvent+fromTrial),...
        '; condition: ', strrep(event{eventIndx}, '_',' '), '; is correct: ',...
        isCorrect{:}, '; RT: ', num2str(rt), ' blinks: ', num2str(bcount)]) % 
%     set(gca,'xlim', [0 3500], 'FontSize',12)
end

if pauseForInspection && plotPlinks
    pause()
end
 

%% find a vector containing boolean values for segment rejection (1 keep the trial)
if ii == size(bIndx,2)
    
a = ceil( blinkPeaks/srate);
a(a == 0) = []; a(a > nrSegments) = [];

seg = zeros(nrSegments,1);     
seg(a,1) = 1;

notblinks = seg < 1;

end
% txt = find(notblinks == 0)

end

if plotPlinks; hold off; end

%% extract some data
if isempty(rt)
    rtPsychoTrials(thisCueEvent+fromTrial) = nan;
else
    rtPsychoTrials(thisCueEvent+fromTrial) = rt;
end

contrastPsychoTrials(thisCueEvent+fromTrial) = contrast;

correctPsychoTrials(thisCueEvent+fromTrial) = str2num(isCorrect{:});

%% prepare fft data   
    if size(elec2plot,2)-1 > 1 
        allSamples=squeeze(mean(allSamples(1:end-1, :))); % mean of channels (not including ocular channel)
    else
        allSamples=squeeze(allSamples(1:end-1, :)); 
    end
    
    % We can now bin into 1 second intervals
    rebinnedData=reshape(allSamples, EEG.srate,trialDur/EEG.srate);
    fftRebinned=fft(rebinnedData); % Perform FFT down time
    % Now the magic happens :) We can average coherently across bins to remove
    % noise
%%
goBack = 0;
%[~, goBack] = findManual(rebinnedData, manualCheck, trialDur/EEG.srate, escapeKey, 2);
% notblinks is a vector where 0 indicates a blink

if sum(notblinks) ~= 0 % If there is at least one segment without a blink
    rebinnedData = rebinnedData(:,find(notblinks));
    
%     blinkIndex = blinkIndex + 1;
%     trialId(thisCueEvent) = thisCueEvent+fromTrial;
%     meanComplexFFTallTrials(:,blinkIndex) = mean(fftRebinned,2); % keep this for all the trials

    meanComplexFFT(:,thisCueEvent) = mean(fftRebinned,2); % nb abs % This is the mean complex FFT for this trial (averaged across bins)
    tIndex(thisCueEvent+fromTrial) = 1;


else  
    fprintf('no good segments found in this trial \n')
end

thisEventBlinks{thisCueEvent,1} = notblinks;
% 

if (thisCueEvent + goBack) > 0
thisCueEvent = thisCueEvent + goBack;
end

%meanComplexFFT(:,thisCueEvent) = mean(fftRebinned,2); % nb abs % This is the mean complex FFT for this trial (averaged across bins)
%thisEventBlinks{thisCueEvent,1} = blinks;   
    
end

% grandAverage{eventIndx}=mean((meanComplexFFT),2); % Average across trials - the abs means we do NOT still keep coherent information
grandAverage{subi, eventIndx} = mean(meanComplexFFT(:,find(tIndex)),2); % Average across trials - the abs means we do NOT still keep coherent information
if ~isempty(meanComplexFFT)   
    allTheSingleTrials{subi, eventIndx} =  abs(meanComplexFFT(1:60,find(tIndex)) ).^2; % Here we lose coherence.
end
allBlinks{subi, eventIndx} = thisEventBlinks;
 

%% SNR

hz = abs(grandAverage{subi, eventIndx}(2:44)).^2;
snrE = zeros(1,size(hz,1));
skipbins =  1; % 1 Hz, hard-coded! (not skipping)
numbins  = 2; %  2 Hz, also hard-coded!

% loop over frequencies and compute SNR
for hzi=numbins+1:length(hz)-numbins-1
    numer = hz(hzi);
    denom = rms( hz([hzi-numbins:hzi-skipbins hzi+skipbins:hzi+numbins]) ); 
    snrE(hzi) = numer./denom;
end  
%%
allSnrE{subi,eventIndx} = snrE;
else
   allSnrE = [];
   allTheSingleTrials = [];
end
end % next event
elecStruct.electrode(eleci).allSnrE = allSnrE;
elecStruct.electrode(eleci).allTheSingleTrials = allTheSingleTrials;
end % next electrode

end % end subject

save([dataDir date '-TFUR-amp-electStruct.mat'], 'elecStruct', '-v7.3');

% save data ressData
% if stimLocked == 1
%     save([dataDir date '-TFUR-SNRs-Occipital-stimLocked.mat'], 'allSnrE', '-v7.3');
%     save([dataDir date '-TFUR-amp-Occipital-stimLocked.mat'], 'allTheSingleTrials', '-v7.3');
%     save([dataDir date '-TFUR-grandAverage-Occipital-stimLocked.mat'], 'grandAverage', '-v7.3');
% else
%     save([dataDir date '-TFUR-SNRs-Occipital-perdLocked-Oz.mat'], 'allSnrE', '-v7.3');
%     save([dataDir date '-TFUR-amp-Occipital-perdLocked-Oz.mat'], 'allTheSingleTrials', '-v7.3');
% end

do.singleElectrodes = 0;
end

%% blinkStructure
% indiDir = 'C:\Users\Richard Naar\Documents\dok\ssvep\Visit to York\EEG data\individual\'; % data after epoching and artefact inspection
dataList = dir([dataDir, '*.mat']); 

% load([dataDir dataList(1).name])
% fprintf('Loading: %s \n', indiList(1).name)

bs = struct();

for subi = 1:size(p.subject,2)
    bs.subject(subi).Id = p.subject(subi).Id;
    for condi = 1:length(p.subject(subi).condition)
        bs.subject(subi).condi(1).blinks = p.subject(subi).condition(4).data.EEG.artefacts;
    end
end

save([dataDir date '-blinkStruct.mat'], 'bs', '-v7.3');

%% Rhythmic Entrainment Source Separation
% only for lower frequency at the moment
while do.ress

event = {'PRED_HIGH' 'PRED_LOW' 'PRED_RND_H' 'PRED_RND_L'};

% frequencies
iLow = [4:4:50];
iHigh = [15:15:50];
imf = [11, 19, 22, 27, 38, 41];
imfTxt = [{'F2 - F1'},{'F1 + F2'},{'(2*F2) - (2*F1)'},{'3*F1 + F2'},{'(F1 + F2)*2'}, {'3*F2 - F1'}];
frex = [iLow, iHigh, imf];
frex = [15, 30, 4, 8, 11];


% parameters for RESS:
srate = 500;
peakwidt  = .5; % FWHM at peak frequency
neighfreq = 1;  % distance of neighboring frequencies away from peak frequency, +/- in Hz
neighwidt = 1;  % FWHM of the neighboring frequencies

% freq1 = 8;
% frex = [8, 30];
% freq2 = 30;
srate = 500;

nfft = ceil( srate/0.1 ); % .1 Hz resolution
hzR    = linspace(0,srate,nfft);

dNr = 3;
matList = dir([dataDir, '*.mat']); 
load([dataDir matList(dNr).name])

fprintf('loading dataset: %s', matList(dNr).name); fprintf('\n') % loading dataset: 24-Nov-2019-TFUR-ress.mat

ressAll = ressAllSegments; % ressAllTrials;% 
clear snrRESS
for subi = 1:size(ressAll.subject,2) % subjects loop
for ressi = 1:size(ressAll.subject(subi).condition,2) % condition loop
    
data_raw = ressAll.subject(subi).condition(ressi).data';
data_raw = data_raw(~cellfun('isempty',data_raw)); % in case some trials got rejected
size(data_raw)

cellDat = {}; data = []; dati = 1;
for ij = 1:size(data_raw,1)
    for uj = 1:size(data_raw{ij},3)
        cellDat(dati,1) = {data_raw{ij}(:,:,uj)};
%         data(:,:, dati) = cellDat{dati,1};
        data = [data, cellDat{dati,1}];

        dati = dati + 1;
    end
end

% data = squeeze(data);
% data = data_raw;
%     
% data = squeeze(ressData(subi,ressi,:))';

%% constantDur determines if all the trials will be stacked together or a constant duration is used
% simulation showed that stacking performs better
% constantDur = 1;
% if constantDur == 1
%     covDat = zeros(size(data{1},1), 750, size(data,2));
%     for tri = 1:size(data,2)
%         covDat(:,:,tri) = data{tri}(:,1:750);
%     end
%     data = covDat;
% % else
%    chan = size(data{1},1); time = size(data{1},2); segments = size(data,1);

   
%    data1 = cell2mat( data ); 
%    data1 = reshape(data1', chan, time, segments);
%     data1 = reshape(data, chan, []);
%
   % end
%%

nbchan = size(data,1);
% data = p.data_stacked ;    
% GED for spatial filter

for frexi = 1:length(frex) % frequencies

[covHi, covLo, covAt] = deal( zeros(nbchan) );
    
% compute covariance matrix for lower neighbor
fdatLo = filterFGx(data,srate,frex(frexi)+neighfreq,neighwidt);
fdatLo = reshape( fdatLo(:,:,:), nbchan,[] );
fdatLo = bsxfun(@minus,fdatLo,mean(fdatLo,2));
covLo  = covLo + (fdatLo*fdatLo')/size(data,2); % pnts?

% compute covariance matrix for upper neighbor
fdatHi = filterFGx(data,srate,frex(frexi)-neighfreq,neighwidt);
fdatHi = reshape( fdatHi(:,:,:), nbchan,[] );
fdatHi = bsxfun(@minus,fdatHi,mean(fdatHi,2));
covHi  = covHi + (fdatHi*fdatHi')/size(data,2);
    
% compute covariance matrix at peak frequency
fdatAt = filterFGx(data,srate,frex(frexi),peakwidt);
fdatAt = reshape( fdatAt(:,:,:), nbchan,[] );
fdatAt = bsxfun(@minus,fdatAt,mean(fdatAt,2));
covAt  = covAt + (fdatAt*fdatAt')/size(data,2);

% perform generalized eigendecomposition. This is the meat & potatos of RESS
[evecs,evals] = eig(covAt,(covHi+covLo)/2);
[~,comp2plot] = max(diag(evals)); % find maximum component

% extract components and force sign
% maps = inv(evecs'); % get maps (this is fine for full-rank matrices)
maps = covAt * evecs / (evecs' * covAt * evecs); % this works either way
[~,idx] = max(abs(maps(:,comp2plot))); % find biggest component
maps = maps * sign(maps(idx,comp2plot)); % force to positive sign

% reconstruct RESS component time series
ress_ts1 = zeros(size(data,2),size(data,3));
for ti=1:size(data,3)
    ress_ts1(:,ti) = evecs(:,comp2plot)'*squeeze(data(:,:,ti));
end    

% ress_ts{subi, ressi} = ress_ts1;
ress_ts = ress_ts1;

%% compute SNR spectrum for RESS

% ressx = mean(abs( fft(ress_ts{subi, ressi}(:,:),nfft,1)/size(data,2) ).^2,2);
ressx = mean(abs( fft(ress_ts,nfft,1)/size(data,2) ).^2,2);
% 
% plot(hzR(1:80/.1), ressx(1:80/.1))
snrR = deal(zeros(size(hzR)));
skipbins =  10; % 10; % 1 Hz, hard-coded!
numbins  = 20+skipbins; % 20+skipbins; %  2 Hz, also hard-coded!


% loop over frequencies and compute SNR
for hzi=numbins+1:length(hzR)-numbins-1
    numer = ressx(hzi);
    denom = mean( ressx([hzi-numbins:hzi-skipbins hzi+skipbins:hzi+numbins]) ); % rms
    snrR(hzi) = numer./denom;
end

snrRESS{subi, ressi, frexi} = snrR;

end % frex
fprintf('doing event: %s', event{ressi}); fprintf('\n')

end
end

save([dataDir date '-allRESSsnr.mat'], 'snrRESS', '-v7.3');
do.ress = 0;
end
%% just to try

implistSNR = dir([dataDir, '*.mat']);  %orderfields(implistSNR, date) 
% [x,idx] = sort(datenum([{implistSNR.date}]), 'descend');
% implistSNR = implistSNR(idx);
dati = 8;
load([dataDir, implistSNR(dati).name]); % load data

fprintf('loading participant: %s \n', implistSNR(dati).name);
%%
srate = 500;
nfft = ceil( srate/0.1 ); % .1 Hz resolution
hzR    = linspace(0,srate,nfft);

meanRESS1 = 0; meanRESS2 = 0; meanRESS3 = 0; meanRESS4 = 0; meanRESS5 = 0;
for subi = 1:23
    freso = 0.1;
    f2show = 120;
    f = 5;
    meanRESS1 = meanRESS1 + snrRESS{subi,1,f}(1:(f2show/freso));
    meanRESS2 = meanRESS2 + snrRESS{subi,2,f}(1:(f2show/freso)); 
    meanRESS3 = meanRESS3 + snrRESS{subi,3,f}(1:(f2show/freso)); 
    meanRESS4 = meanRESS4 + snrRESS{subi,4,f}(1:(f2show/freso)); 
%     meanRESS5 = meanRESS5 + snrRESS{subi,1,2}(1:(f2show/freso)); 
    
end
hz2plot = hzR(1:(f2show/freso));

stem(hz2plot, meanRESS1/23, 'k')
hold
stem(hz2plot, meanRESS2/23, 'r')
stem(hz2plot, meanRESS3/23, 'g')
stem(hz2plot, meanRESS4/23, 'b')
%plot(hz2plot, meanRESS5/23, 'y')
hold off
% legend({'15';'30';'4'; '8'; '11'})
legend({'high';'low';'high50'; 'low50'})

%% 11 19 22 29 44 (colFun((15^2)+(4^2), 1) and colFun((15^2)-(4^2), 1);)
sub = 1;
 freso = 0.1;
% for subi = 1:size(snrRESS,1)
    for fi = 1:size(snrRESS,3)
        f = fi; sub = subi;
        f2show = 40 ;
        plot(hzR(1:(f2show/freso)), snrRESS{sub,1,f}(1:(f2show/freso)), 'k')
        hold
        plot(hzR(1:(f2show/freso)), snrRESS{sub,2,f}(1:(f2show/freso)), 'r')
        plot(hzR(1:(f2show/freso)), snrRESS{sub,3,f}(1:(f2show/freso)),'--g')
        plot(hzR(1:(f2show/freso)), snrRESS{sub,4,f}(1:(f2show/freso)),'--b')
        title(['sub:', num2str(subi) ,'; frquency:' ,num2str(frex(f)), ' Hz'])
        legend({'High';'Low';'High RND'; 'Low RND'})
        hold
        pause
    end
% end
%% put the data into table 1 (allSnrE)

while do.saveCsv == 1
matList = dir([dataDir, '*.mat']); 
dati = 16;
load([dataDir matList(dati).name])

fprintf('loading dataset: %s', matList(dati).name); fprintf('\n')

%allSnrE = snrRESS;

iLow = [4:4:40];
iHigh = [15:15:40];
imf = [11, 19, 22, 27, 38];
% imfTxt = [{'F2 - F1'},{'F1 + F2'},{'(2*F2) - (2*F1)'},{'3*F1 + F2'},{'(F1 + F2)*2'}, {'3*F2 - F1'}];
allFrex = [iLow, iHigh, imf];

% FOIsLow = zeros(size(subtraction)); FOIsLow(iLow) = subtraction(iLow);

header = [{'SubId', 'cond',...
    'low4', 'low8', 'low12', 'low16', 'low20', 'low24','low28', 'low32', 'low36', 'low40',...
    'high15', 'high30', ...
    'F11', 'F19', 'F22','F27', 'F38'}];

data = {};
data(1, :) = header;
subDat = [];
for rind = 1:size(allSnrE,1)
%   for win = 1:length(allFrex)
%       currentWin = allFrex(win)-0.5:allFrex(win)+0.5
%       high = allSnrE{rind,1}(currentWin); % PRED_HIGH
       
       
%V1       
      high = allSnrE{rind,1}(allFrex); % PRED_HIGH
      low =  allSnrE{rind,2}(allFrex); % PRED_LOW
      high_rnd = allSnrE{rind,3}(allFrex); % PRED_RND_H
      low_rnd = allSnrE{rind,4}(allFrex); % PRED_RND_L


%       high_rnd = allSnrE{rind,3}(allFrex); % PRED_RND_H
%       low_rnd = allSnrE{rind,4}(allFrex); % PRED_RND_L
%       rnd_mean = (high_rnd + low_rnd)/2;
%       
%       high = allSnrE{rind,1}(allFrex) - rnd_mean; % PRED_HIGH
%       low =  allSnrE{rind,2}(allFrex) - rnd_mean; % PRED_LOW

      
      id = num2str(rind);
%   end
  if rind < 2
      top = 1; bottom = 4; %V1
%       top = 1; bottom = 2;
  else
      top = top+3 ; bottom = bottom+3; %V1
%       top = top+1 ; bottom = bottom+1; 
  end
%V1  
  data(rind+top:rind+bottom, :) =  [{id}, {'Cued high'}, num2cell(high); {id}, {'Cued low'}, num2cell(low); {id},...
      {'Non-cued high'}, num2cell(high_rnd); {id}, {'Non-cued low'}, num2cell(low_rnd)];
%   data(rind+top:rind+bottom, :) =  [{id}, {'Cued high'}, num2cell(high); {id}, {'Cued low'}, num2cell(low)];
   
end

addpath('C:\Users\Richard Naar\Documents\Matlab toolboxes\letswave6-master\external')

cd(dataDir)
cell2csv([date '-TFUR-SNRs-Pz-Poz-Oz.csv'], data)

do.saveCsv = 0;
end

%% put the data into table 2 (snrRESS)

while do.saveCsvRess == 1
matList = dir([dataDir, '*.mat']); 
dati = 9;
load([dataDir matList(dati).name])

fprintf('loading dataset: %s', matList(dati).name); fprintf('\n')

iLow = [4:4:50];
iHigh = [15:15:50];
imf = [11, 19, 22, 27, 38, 41];
imfTxt = [{'F2 - F1'},{'F1 + F2'},{'(2*F2) - (2*F1)'},{'3*F1 + F2'},{'(F1 + F2)*2'}, {'3*F2 - F1'}];
allFrex = [iLow, iHigh, imf]/0.1;
allFrex = [15, 30, 4, 8, 11]/0.1;
snrDat = snrRESS;
% FOIsLow = zeros(size(subtraction)); FOIsLow(iLow) = subtraction(iLow);

header = [{'SubId', 'cond',...
         'low4', 'low8',...
         'high15', 'high30',...
         'F11'}];
data = {};
data(1, :) = header;
subDat = [];
for rind = 1:size(snrDat ,1)
[high, low, high_rnd, low_rnd] = deal(zeros(1, size(snrDat ,3)));    
for frexi = 1:size(snrDat ,3)
    
  high(1,frexi) = max(snrDat{rind,1,frexi}(allFrex(frexi)-2:allFrex(frexi)+2)); % PRED_HIGH
  low(1,frexi) =  max(snrDat{rind,2,frexi}(allFrex(frexi)-2:allFrex(frexi)+2)); % PRED_LOW
  high_rnd(1,frexi) = max(snrDat{rind,3,frexi}(allFrex(frexi)-2:allFrex(frexi)+2)); % PRED_RND_H
  low_rnd(1,frexi) = max(snrDat{rind,4,frexi}(allFrex(frexi)-2:allFrex(frexi)+2)); % PRED_RND_L
  id = num2str(rind);
  
end


top = size(data,1)+1; bottom = size(data,1)+4;

data(top:bottom, :) =  [{id}, {'Cued high'}, num2cell(high); {id}, {'Cued low'}, num2cell(low); {id},...
      {'Non-cued high'}, num2cell(high_rnd); {id}, {'Non-cued low'}, num2cell(low_rnd)];


end

addpath('C:\Users\Richard Naar\Documents\Matlab toolboxes\letswave6-master\external')

cd(dataDir)
cell2csv([date '-TFUR-SNRs-RESS-occipital-pred.csv'], data)

do.saveCsvRess = 0;
end

%% put the data into table

while do.saveCsvSingle == 1
matList = dir([dataDir, '*.mat']); 
load([dataDir matList(4).name])

fprintf('loading dataset: %s', matList(4).name); fprintf('\n')

iLow = [4:4:50];
iHigh = [15:15:50];
imf = [11, 19, 22, 27, 38, 41];
imfTxt = [{'F2 - F1'},{'F1 + F2'},{'(2*F2) - (2*F1)'},{'3*F1 + F2'},{'(F1 + F2)*2'}, {'3*F2 - F1'}];
allFrex = [iLow, iHigh, imf];

% addpath('C:\Users\Richard Naar\Documents\fNIRS\Kersteni uurimistöö\data\analysisR\fnirs')

% FOIsLow = zeros(size(subtraction)); FOIsLow(iLow) = subtraction(iLow);

header = [{'SubId', 'cond',...
    'low4', 'low8', 'low12', 'low16', 'low20', 'low24','low28', 'low32', 'low36', 'low40', 'low44', 'low48',...
    'high15', 'high30', 'high45',...
    'F11', 'F19', 'F22','F27', 'F38','F41', 'trial'}];
data = {};
data(1, :) = header;
subDat = [];


for rind = 1:size(allTheSingleTrials,1)
    
%% SNR
nTrials = zeros(1,4);
clear allTrialsSnrE
for eventi = 1:4
nTrials(eventi) = size(allTheSingleTrials{rind,eventi},2);        
for triali = 1:nTrials(eventi)    
hz = allTheSingleTrials{rind, eventi}(2:end, triali);
snrE = zeros(1,size(hz,1));
skipbins =  0; % 1 Hz, hard-coded! (not skipping)
numbins  = 1; %  2 Hz, also hard-coded!

% loop over frequencies and compute SNR
for hzi=numbins+1:length(hz)-numbins-1
    numer = hz(hzi);
    denom = rms( hz([hzi-numbins:hzi-skipbins hzi+skipbins:hzi+numbins]) ); 
    snrE(hzi) = numer./denom;
end  
%%
% %     subplot(2,2,eventIndx); hold on
% %     title(strrep(event{eventIndx}, '_',' '))

%    bar(abs(grandAverage(2:60)));

%     bar(abs(grandAverage{eventIndx}(2:60)));
allTrialsSnrE{1,eventi}(:, triali) = snrE;
end
end
%% 
    
  high = allTrialsSnrE{1}(allFrex,:); % PRED_HIGH
  low =  allTrialsSnrE{2}(allFrex,:); % PRED_LOW
  high_rnd = allTrialsSnrE{3}(allFrex,:); % PRED_RND_H
  low_rnd = allTrialsSnrE{4}(allFrex,:); % PRED_RND_L

%   nTrials = size(allTheSingleTrials{rind,1},2);
%     
%   high = allTheSingleTrials{rind,1}(allFrex,:); % PRED_HIGH
%   low =  allTheSingleTrials{rind,2}(allFrex,:); % PRED_LOW
%   high_rnd = allTheSingleTrials{rind,3}(allFrex,:); % PRED_RND_H
%   low_rnd = allTheSingleTrials{rind,4}(allFrex,:); % PRED_RND_L
  
% for trial = 1:size(low_rnd,2)
  
  id = num2str(rind);
  
  if rind < 2
      top = 1; bottom = sum([size(high,2), size(low,2), size(high_rnd,2), size(low_rnd,2)]);
      addVal = 1;
  else
      addVal = 0;
      currentSub = sum([size(high,2), size(low,2), size(high_rnd,2), size(low_rnd,2)]);
      top = size(data,1)+1 ; bottom = (top+currentSub)-1;
  end
  
%   data(rind+addVal1:rind+addVal2, :) =  [repmat([{1}, {'Cued high'}],nTrials(1) ,1), num2cell(high)', num2cell(1:nTrials(1))'];
%     
  
  data(top+addVal:bottom+addVal, :) =  [repmat([{id}, {'Cued high'}],nTrials(1) ,1), num2cell(high)', num2cell(1:nTrials(1))'; repmat([{id}, {'Cued low'}],nTrials(2) ,1), num2cell(low)', num2cell(1:nTrials(2))'; ...
      repmat([{id}, {'Non-cued high'}],nTrials(3) ,1), num2cell(high_rnd)', num2cell(1:nTrials(3))'; repmat([{id}, {'Non-cued low'}],nTrials(4) ,1), num2cell(low_rnd)', num2cell(1:nTrials(4))'];

  
% end
end

addpath('C:\Users\Richard Naar\Documents\Matlab toolboxes\letswave6-master\external')

cd(dataDir)
cell2csv([date '-TFUR-SingleTrials-SNR-PRED.csv'], data)

do.saveCsvSingle = 0;
end

%% plotting

while do.plotting == 1
% 11, 19, 22, 27,38     
% low = 4;
% high = 15;
% imf = low+high
% imf2 = high - low
% imf*2 - 11

% nF + m f
%% plot frex

event = {'PRED_HIGH' 'PRED_LOW' 'PRED_RND_H' 'PRED_RND_L'}; 
% event = {'STIM_HIGH' 'STIM_LOW' 'STIM_RND_H' 'STIM_RND_L'};
% event = {'CUE_HIGH' 'CUE_LOW' 'CUE_RND_H' 'CUE_RND_L'};

implistSNR = dir([dataDir, '*.mat']);  %orderfields(implistSNR, date) 
% [x,idx] = sort(datenum([{implistSNR.date}]), 'descend');
% implistSNR = implistSNR(idx);

%% load data and trim frex
dati = 5;%STIM: 28; 5; 6; %CUE: 15;% 16; %17; TFUR-grandAverage-OC-stimLocked

load([dataDir, implistSNR(dati).name]); % load data
fprintf('loading participant: %s \n', implistSNR(dati).name);

allSnrE = grandAverage;

%     nFrex = size([allSnrE{1,1}],2); %nFrex = size([allSnrE{1,1}],1);
nFrex = 43; %nFrex = size([allSnrE{1,1}],1);

%     hz = abs(grandAverage{subi, eventIndx}(2:44)).^2;

datStandard = [];
for condi = 1:4% size(allSnrE,2)
    datn = squeeze(cat(3,allSnrE{:,condi}));
    datn = datn(2:nFrex+1, :); % excluding 0 Hz
    datStandard = [datStandard; datn];
end

%  currently not standardizing 
% meanCentered = bsxfun(@minus,datStandard, mean(datStandard));
% modifZdat = bsxfun(@rdivide, meanCentered, mean(datStandard));
% medianCentered = .6745.*bsxfun(@minus,datStandard,median(datStandard));
% modifZdat = bsxfun(@rdivide, medianCentered, median(abs(medianCentered)));
modifZdat = datStandard;


%%
for frex = 1:2
comp1 = frex;  
% comp2 = 3;
figN = 1;
% figure(figN)
foiLow = [4:4:40];  %foiLow = foiLow+1;
foiHigh = [15:15:40];%  foiHigh = foiHigh+1;
% high
figure(frex)
subplot(3,1,1); hold on
% bar(grandAverage{comp1}(2:60));
% bar(mean( squeeze(cat(3,allSnrE{:,comp1})), 2));

if comp1 == 1
    predDat = abs(mean(modifZdat(1:nFrex,:),2)).^2; % high 
else
    predDat = abs(mean(modifZdat(nFrex+1:nFrex*2,:),2)).^2; % low
end
    
bar(predDat)

imf = [11, 19, 22, 27, 38, 41];
% imfTxt = ['F2-F1','F1+F2','2*F2-2*F1','3*F1+F2','(F1+F2)*2'];
imfTxt = [{'F2 - F1'},{'F1 + F2'},{'(2*F2) - (2*F1)'},{'3*F1 + F2'},{'(F1 + F2)*2'}, {'3*F2 - F1'}];
% maxY = max(mean( squeeze(cat(3,allSnrE{:,comp1})), 2))+1;
maxY = max(predDat)+100; minY = min(predDat);
% text((imf)-1, zeros(1,length(imf(1:length(imf))))+maxY, imfTxt(1:length(imf))) % *mean(EEG.data(31,:))
% set(gca,'ylim',[min(mean( squeeze(cat(3,allSnrE{:,comp1})), 2)) maxY+2], 'FontSize',12)
% set(gca,'ylim',[minY maxY], 'FontSize',12)


%% plot the colors
% FOIsLow = zeros(size(grandAverage{1})); FOIsLow(foiLow) = grandAverage{1}(foiLow);
% bar(FOIsLow(2:60), 'r')
% % color high
% FOIsHigh = zeros(size(grandAverage{1})); FOIsHigh(foiHigh) = grandAverage{1}(foiHigh);
%% bar(FOIsHigh(2:60), 'g')
title(strrep(event{comp1}, '_',' '))
%% high - random
% subtraction = grandAverage{comp1}(1:60)-(grandAverage{3}(1:60)+grandAverage{4}(1:60))/2;
% subtraction_lrnd = mean( squeeze(cat(3,allSnrE{:,comp1})), 2)-(mean( squeeze(cat(3,allSnrE{:,3})), 2)+mean( squeeze(cat(3,allSnrE{:,4})), 2))/2;

highrnd = mean(modifZdat(nFrex*2+1:nFrex*3,:),2);
lowrnd = mean(modifZdat(nFrex*3+1:nFrex*4,:),2);

subtraction_lrnd = predDat - (lowrnd+highrnd)/2;



subplot(3,1,2); hold on
bar(subtraction_lrnd);
maxY = max(subtraction_lrnd)+10; minY = min(subtraction_lrnd)-10;
% text((imf)-1, zeros(1,length(imf(1:length(imf))))+maxY, imfTxt(1:length(imf)))
% set(gca,'ylim',[min(subtraction)-2 maxY+2], 'FontSize',12)
%% plot the colors
% color low

FOIsLow = zeros(size(subtraction_lrnd)); FOIsLow(foiLow) = subtraction_lrnd(foiLow);
bar(FOIsLow, 'm')

% color high
FOIsHigh = zeros(size(subtraction_lrnd)); FOIsHigh(foiHigh) = subtraction_lrnd(foiHigh);
bar(FOIsHigh, 'k')
title('certain - uncertain condition')
legend({'Other';'Low';'High'})
% set(gca,'ylim',[minY maxY], 'FontSize',12)

%% random high
subplot(3,1,3);
avgRnd = (lowrnd+highrnd)/2;
avgRand = subtraction_lrnd;
bar(avgRnd);
maxY = max(avgRnd)+1;
% text((imf)-1, zeros(1,length(imf(1:length(imf))))+maxY, imfTxt(1:length(imf)))
% set(gca,'ylim',[min(avgRnd) maxY+2], 'FontSize',12)
%  title(strrep(event{comp2}, '_',' '))
title('RANDOM')
pause()
end
%% both on the same
% for subi = 1:23
% collatz = [30, 15, 46, 23, 70, 35, 106, 53, 160, 80, 40, 20, 10, 5, 16, 8, 4, 2, 1]
figure(1)
leg = 1;
comp1 = 1;  
maxF = 40;

% color
colHigh = [0 0 .7];
colLow = [.7 0 0];
imfCol = [0 .7 0];
colPins = [0 0 0];

figN = 1;
% figure(figN)
foiLow = [4:4:maxF];  
foiHigh = [15:15:maxF];

% higher frequency
% subplot(2,1,1); 
hold on
% highDat = squeeze(cat(3,allSnrE{:,comp1})); highDat = highDat(1:maxF,:);
% randoDat = (squeeze(cat(3,allSnrE{:,3}))+squeeze(cat(3,allSnrE{:,4}))) / 2; randoDat = randoDat(1:maxF,:);

% lowDat = abs(mean(modifZdat(nFrex+1:nFrex*2,:),2)).^2;
% highDat = abs(mean(modifZdat(1:nFrex,:),2)).^2;

% bsxfun(@minus, modifZdat(1:nFrex,:), randoDat)

highDat = abs(modifZdat(1:nFrex,:)).^2; % highDat = highDat/1000; % convert to mV
highDatMean = mean(highDat,2);
highDat_numer = bsxfun(@minus, highDat, min(highDatMean));
highDat = bsxfun(@rdivide, highDat_numer, max(highDatMean) - min(highDatMean));

lowDat = abs(modifZdat(nFrex+1:nFrex*2,:)).^2; % lowDat = lowDat/1000; % convert to mV
lowDatMean = mean(lowDat,2);
lowDat_numer = bsxfun(@minus, lowDat, min(lowDatMean ));
lowDat = bsxfun(@rdivide, lowDat_numer, max(lowDatMean ) - min(lowDatMean ));

% highrnd = abs(mean(modifZdat(nFrex*2+1:nFrex*3,:),2)).^2;
% lowrnd = abs(mean(modifZdat(nFrex*3+1:nFrex*4,:),2)).^2;

% randoDat = (lowrnd+highrnd)/2;
% randoDat = squeeze(cat(3,allSnrE{:,comp1+1})); randoDat = randoDat(1:maxF,:);

% subtraction_hl = (mean( highDat, 2) - mean(lowDat, 2))./(mean( highDat, 2) + mean(lowDat, 2));
subtraction_hl = mean( highDat, 2) - mean(lowDat, 2);

% median
% % highDat = median(modifZdat(1:nFrex,:),2);
% % highrnd = median(modifZdat(nFrex*2+1:nFrex*3,:),2);
% % lowrnd = median(modifZdat(nFrex*3+1:nFrex*4,:),2);
% % 
% % 
% % randoDat = (lowrnd+highrnd)/2;
% % % randoDat = squeeze(cat(3,allSnrE{:,comp1+1})); randoDat = randoDat(1:maxF,:);
% % 
% % subtraction_hrnd = median( highDat, 2) - median(randoDat, 2);



bar(subtraction_hl, 'FaceColor', colPins);
imf = [11, 19, 22, 27, 38, 41];
imf = imf(find(imf < maxF));

% imfTxt = [{'F2 - F1'},{'F1 + F2'},{'(2*F2) - (2*F1)'},{'3*F1 + F2'},{'(F1 + F2)*2'}, {'3*F2 - F1'}];
% maxY = max(subtraction_hl)+1;
% text((imf)-1, zeros(1,length(imf(1:length(imf))))+maxY, imfTxt(1:length(imf))) % *mean(EEG.data(31,:))

%% plot the colors
% color low - high

FOIsLow = zeros(size(subtraction_hl)); FOIsLow(foiLow) = subtraction_hl(foiLow);
bar(FOIsLow, 'FaceColor', colLow)

% color high
FOIsHigh = zeros(size(subtraction_hl)); FOIsHigh(foiHigh) = subtraction_hl(foiHigh);
bar(FOIsHigh, 'FaceColor', colHigh)
% title([strrep(event{comp1+2}, '_',' '), ' - ',strrep(event{comp1+3}, '_',' ')])

% color imf 
FOIsImf = zeros(size(subtraction_hl)); FOIsImf(imf) = subtraction_hl(imf);
bar(FOIsImf, 'FaceColor', imfCol)

if leg; legend({'Other';'Low';'High';'Imf'}); end

% bar(FOIsHigh(2:60), 'g')
% title([strrep(event{comp1}, '_',' '), ' - ',strrep(event{comp1+1}, '_',' ')])
% title('High-Low')

% 
% err
% mDat = mean( highDat, 2) - mean(randoDat, 2); 
% 
% % highDat = mean(modifZdat(1:nFrex,:),2);
% % highrnd = mean(modifZdat(nFrex*2+1:nFrex*3,:),2);
% % lowrnd = mean(modifZdat(nFrex*3+1:nFrex*4,:),2);
% % 
% % randoDat = (lowrnd+highrnd)/2;
% 
% scaledStd = (std(bsxfun(@minus, modifZdat(1:nFrex,:), randoDat),0,2)./ sqrt(22) )  * 1.96;
% 
% % scaledStd = ( std( highDat - randoDat  ,0,2) ./ sqrt(22) )  * 1.96;
% errHigh = mDat + scaledStd;
% errLow = mDat - scaledStd;
% %
% k = 0.2; ls = 1.5;
% for pli = 1:length(mDat)
%     plot([pli-k pli+k] ,[errHigh(pli) errHigh(pli)], '-', 'LineWidth',ls,'Color', 'k') % upper vertical
%     plot([pli pli] ,[errHigh(pli) errLow(pli)], '-','LineWidth',ls, 'Color', 'k') % horizontal
%     plot([pli-k pli+k] ,[errLow(pli) errLow(pli)], '-','LineWidth',ls, 'Color', 'k') % lower vertical
% end
% set(gca,'ylim',[floor(min(errLow))-5 ceil(max(errHigh))+5], 'FontSize',12)

% set(gca,'ylim',[-5 30])

%% 
subtr = highDat - lowDat;
% subtr2 = (highDat - lowDat)./(highDat + lowDat);


N = 23;
err = std(subtr')./sqrt(N);


errorb(1:43,subtraction_hl,err,'top')


% % 
%%
% set(gca,'ylim',[-5 30], 'FontSize',12)
%%%set(gca,'xlim',[1 40], 'FontSize',12)
ylabel('Amplitude (mV)'); xlabel('Frequency')  % Modified z-scores
set(findall(gcf,'-property','FontSize'),'FontSize',14)
set(findall(gcf,'-property','fontname'),'fontname', 'Arial')


%% sig
% 
labdat = subtraction_hl+err';
sigf = [30];
% text(sigf-0.15, labdat(sigf)+1, '*', 'fontsize',12)    
text(sigf-1, labdat(sigf)+1, 'p = .052', 'fontsize',12)    


hold off


%% small plot


fs = 29:31;
% fs = 7:9;
figure(2)
hold on
bar(fs, subtraction_hl(fs), 'FaceColor', colPins);
bar(fs, FOIsImf(fs), 'FaceColor', imfCol)
bar(fs, FOIsHigh(fs), 'FaceColor', colHigh)
bar(fs, FOIsLow(fs), 'FaceColor', colLow)

errorb(fs, subtraction_hl(fs),err(fs),'top')

ylabel('Amplitude (mV)'); xlabel('Frequency')  % Modified z-scores
title('L')
set(findall(gcf,'-property','FontSize'),'FontSize',14)
set(findall(gcf,'-property','fontname'),'fontname', 'Arial')

% % sig
% labdat = subtraction_hl(fs)+err(fs)';
% sigf = [2];
% text(fs(sigf)-0.1, labdat(sigf)+0.2, '*', 'fontsize',18)    
% text(fs(sigf)-0.3, labdat(sigf)+0.5, 'p = 0.052', 'fontsize',10)    


x0=500;
y0=200;
width=300;
height=250;
set(gcf,'position',[x0,y0,width,height])
% cue
% set(gca,'ylim',[-1.5 3])
% set(gca,'ylim',[-3 6])

% stim
% set(gca,'ylim',[-2.5 5])
set(gca,'ylim',[-0.75 2])
%% AVG RND
% lowDat = squeeze(cat(3,allSnrE{:,2})); lowDat = lowDat(1:maxF,:);
% lowDat = abs(mean(modifZdat(nFrex+1:nFrex*2,:),2)).^2;

% lowrnd = abs(mean(modifZdat(nFrex*3+1:nFrex*4,:),2)).^2;
% highrnd = abs(mean(modifZdat(nFrex*2+1:nFrex*3,:),2)).^2;

lowrnd = modifZdat(nFrex*3+1:nFrex*4,:);
highrnd = modifZdat(nFrex*2+1:nFrex*3,:);

rnd = abs((lowrnd + highrnd)/2).^2;
rnd = rnd/1000; % convert to mV

% lowrnd = abs(modifZdat(nFrex*3+1:nFrex*4,:)).^2;
% highrnd = abs(modifZdat(nFrex*2+1:nFrex*3,:)).^2;

rnds = mean(rnd,2);%  mean( lowDat, 2)- mean( randoDat, 2); % these names are wrong

%median
% % lowDat = median(modifZdat(nFrex+1:nFrex*2,:),2);
% % 
% % subtraction_lrnd = median( lowDat, 2)- median( randoDat, 2);

% subplot(2,1,2); 
hold on
bar(rnds, 'FaceColor', colPins);

maxY = max(rnds)+1;
% text((imf)-1, zeros(1,length(imf(1:length(imf))))+maxY, imfTxt(1:length(imf)))
% set(gca,'ylim',[min(subtraction)-2 maxY+2], 'FontSize',12)
%% plot the colors
% color low

FOIsLow = zeros(size(rnds)); FOIsLow(foiLow) = rnds(foiLow);

bar(FOIsLow, 'FaceColor', colLow)


% color high
FOIsHigh = zeros(size(rnds)); FOIsHigh(foiHigh) = rnds(foiHigh);
bar(FOIsHigh, 'FaceColor', colHigh)

% imfTxt = [{'F2 - F1'},{'F1 + F2'},{'(2*F2) - (2*F1)'},{'3*F1 + F2'},{'(F1 + F2)*2'}, {'3*F2 - F1'}];
% maxY = max(subtraction)+1;
% text((imf)-1, zeros(1,length(imf(1:length(imf))))+maxY, imfTxt(1:length(imf))) % *mean(EEG.data(31,:))

% color imf 
FOIsImf = zeros(size(rnds)); FOIsImf(imf) = rnds(imf);
bar(FOIsImf, 'FaceColor', imfCol)

% legend({'Other';'Low';'High';'Imf'})

% err
% mDat = mean( lowDat, 2) - mean(randoDat, 2); 
% 
% scaledStd = (std(bsxfun(@minus, modifZdat(nFrex+1:nFrex*2,:), randoDat),0,2)./ sqrt(22) )  * 1.96;
% 
% % scaledStd = ( std( lowDat - randoDat  ,0,2) ./ sqrt(22) )  * 1.96;
% errHigh = mDat + scaledStd;
% errLow = mDat - scaledStd;
% 
% k = 0.2; ls = 1.5;
% for pli = 1:length(mDat)
%     plot([pli-k pli+k] ,[errHigh(pli) errHigh(pli)], '-', 'LineWidth',ls,'Color', 'k') % upper vertical
%     plot([pli pli] ,[errHigh(pli) errLow(pli)], '-','LineWidth',ls, 'Color', 'k') % horizontal
%     plot([pli-k pli+k] ,[errLow(pli) errLow(pli)], '-','LineWidth',ls, 'Color', 'k') % lower vertical
% end
% set(gca,'ylim',[floor(min(errLow))-5 ceil(max(errHigh))+5], 'FontSize',12)
% 
% 
% %%%set(gca,'xlim',[1 40], 'FontSize',12)
% % set(gca,'ylim',[-30 15], 'FontSize',12)
% 
ylabel('Amplitude (mV)', 'FontSize',14); xlabel('Frequency')  
% title('Attending low')
% title([strrep(event{comp1+2}, '_',' '), ' - ',strrep(event{comp1+3}, '_',' ')])
set(gca,'ylim',[min(rnds)-(min(rnds)*0.1) max(rnds)+(max(rnds)*0.5)])

set(findall(gcf,'-property','FontSize'),'FontSize',14)
set(findall(gcf,'-property','fontname'),'fontname', 'Arial')
% set(gca, 'fontname', 'Arial') % get(gca,'fontname') listfonts

legend({'Other';'Low';'High';'Imf'})

%% err



N = 23;
% err = 1.96 * (std(rnd')./sqrt(N));
err = std(rnd')./sqrt(N);

errorb(1:43,mean(rnd,2),err, 'top')


%%
% 
% juku = mean(rnd,2)+err';
% sigf = [8, 11, 19, 30]
% text(sigf-0.15, juku(sigf)+1000, '*', 'fontsize',18)    

hold off

%% small plot


fs = 29:31;
% fs = 7:9;
figure(2)
hold on
bar(fs, rnds(fs), 'FaceColor', colPins);
bar(fs, FOIsImf(fs), 'FaceColor', imfCol)
bar(fs, FOIsHigh(fs), 'FaceColor', colHigh)
bar(fs, FOIsLow(fs), 'FaceColor', colLow)

errorb(fs, rnds(fs),err(fs),'top')

ylabel('Amplitude (mV)'); xlabel('Frequency')  % Modified z-scores
title('R')
set(findall(gcf,'-property','FontSize'),'FontSize',14)
set(findall(gcf,'-property','fontname'),'fontname', 'Arial')

% % sig
% labdat = rnds(fs)+err(fs)';
% sigf = [2];
% text(fs(sigf)-0.1, labdat(sigf)+1, '*', 'fontsize',18)    
% text(fs(sigf)-0.3, labdat(sigf)+0.5, 'p = 0.052', 'fontsize',10)    


x0=500;
y0=200;
width=300;
height=250;
set(gcf,'position',[x0,y0,width,height])
% cue
% set(gca,'ylim',[0 40])
% set(gca,'ylim',[0 7])

%%


% pause()
end
%%

%prelim stat

% juku = rnd;
% juku = bsxfun(@minus, abs(modifZdat(1:nFrex,:)).^2, abs(modifZdat(nFrex+1:nFrex*2,:)).^2);
juku = bsxfun(@minus, highDat, lowDat);
% juku = subtr2;


% hist(modifZdat(8, :))
% 
% raw
% highrndraw = mean(datStandard(nFrex*2+1:nFrex*3,:),2);
% lowrndraw = mean(datStandard(nFrex*3+1:nFrex*4,:),2);
% randoDatRaw = (lowrndraw+highrndraw)/2;
% juku = bsxfun(@minus, datStandard(nFrex+1:nFrex*2,:), randoDatRaw);



% highDat(30)*sqrt(23)/std(juku(30, :))

frex2comp = 22;
% data1 = (juku(frex2comp-1,:)+juku(frex2comp+1,:))/2;
data1 = zeros(23,1);% juku(7,:);
data2 = juku(frex2comp,:);
% hist(data2)
% [H, pValue, SWstatistic] = swtest(data2);
% pValue

figure(2), clf, hold on

colors = 'kr';
N=max(size(data1));
for i=1:N
    plot([data1(i) data2(i)],[i i],colors((data1(i)<data2(i))+1),'HandleVisibility','off')
end



plot(data1,1:N,'ks','markerfacecolor','k','markersize',10)
plot(data2,1:N,'ro','markerfacecolor','r','markersize',10)

% set(gca,'ylim',[0 N+1],'xlim',[-.5 max([data2])+.5],'xtick',0:5)
ylabel('Data index'), xlabel('Data value')
% grid minor
legend({'data1';'data2'})


[H, pValue, SWstatistic]  = swtest(data2);

if pValue < 0.05
    [p,h,stats] = signrank(data2);
    title([ 'Wilcoxon z=' num2str(stats.zval) ', p=' num2str(p) ])
else
    [h,p,c,stats] = ttest(data2);
    title([ 't=' num2str(stats.tstat) ', p=' num2str(p) ])
end
% end

p
stats

r = stats.zval/sqrt(N);
fprintf('The effect size r is: %f \n', r)
%% PALAMEDES
dataDir = 'C:\Users\Richard Naar\Documents\dok\ssvep\Visit to York\EEG data\';  % directory for datafiles
addpath('C:\Users\Richard Naar\Documents\MATLAB\palamedes1_10_3\Palamedes')

%% load SNR data 

event = {'PRED_HIGH' 'PRED_LOW' 'PRED_RND_H' 'PRED_RND_L'}; 
% event = {'STIM_HIGH' 'STIM_LOW' 'STIM_RND_H' 'STIM_RND_L'};
% event = {'CUE_HIGH' 'CUE_LOW' 'CUE_RND_H' 'CUE_RND_L'};

implistSNR = dir([dataDir, '*.mat']);  %orderfields(implistSNR, date) 
% [x,idx] = sort(datenum([{implistSNR.date}]), 'descend');
% implistSNR = implistSNR(idx);

load([dataDir, implistSNR(8).name]); % load data

fprintf('loading participant: %s \n', implistSNR(8).name);
%% load the blinkStructure

% load([dataDir, implistSNR(15).name]); % load data
% 
% fprintf('loading participant: %s \n', implistSNR(15).name);
% 
% trialsLeft = zeros( size(bs.subject), length(event));
% 
% kuku = bs.subject(1).condi.blinks{1,1};
% kuku = cell(kuku);
% cond1 = cell2mat(kuku);
%% load in the level averages
cd('C:\Users\Richard Naar\Documents\dok\ssvep\Visit to York\EEG data')

psyList = dir([dataDir, '*.txt']);  %orderfields(implistSNR, date) 
%load([dataDir, psyList(14).name]); % load data

% meanContrast = importfile(psyList(5).name, 2, 26);
meanContrast = importMyFile(psyList(3).name); % 2 5981
% meanContrast.Properties.VariableNames = {'rowi' 'participant' 'contrast' 'trials_2label' 'NumPos' 'N'};
% meanContrast.Properties.VariableNames = {'rowi' 'contrast' 'trials_2label' 'NumPos' 'N'};


meanContrast.Properties.VariableNames = {'rowi' 'contrast' 'trials_2label' 'participant' 'trials_2thisRepN' 'trials_2intensity' 'NumPos' 'N'};
[C, IA, IC] = unique(meanContrast.trials_2label);
% conds = {'Cued High' 'Cued Low' 'Non-Cued (high)' 'Non-Cued (low)'};
conds = {meanContrast.trials_2label{IA(2)}, meanContrast.trials_2label{IA(1)}, ...
    meanContrast.trials_2label{IA(3)}, meanContrast.trials_2label{IA(4)}};
%% now data for each subject
plotting = 1;
pictureBase = 'C:\Users\Richard Naar\OneDrive - Tartu Ülikool\Tempral frequency manuscript\Muu\Pildid\individual';

frex = [8,30,11];
cols = {'r' 'g' 'b' 'k'};
colsp = {'.r' '.g' '.b' '.k'};
% 
%for subi = 1:23;
subi = 1;
figure(subi)
hold on
% for plotting = 1:2
for condi = 1:4;
% condi = 3;
% frexi = 1;

%% extract data
%EEG
%dataSnr(subi, condi) = allSnrE{subi,condi}(frex(frexi)); 
% condi = 1;
% psycoph
% condi = 1;
subCondTable = meanContrast(strcmp(meanContrast.trials_2label, conds{condi}) , :); %& (meanContrast.participant == subi )
%subCondTable = meanContrast(strcmp(meanContrast.trials_2label, conds{condi}) & (meanContrast.participant == subi ) , :); %

OutOfNum = subCondTable.N';
% StimLevels = subCondTable.contrast'; % subCondTable.contrast/100;
StimLevels = subCondTable.trials_2intensity'; 
NumPos = subCondTable.NumPos';

%% Fitting psychometric function useing a maximum likelihood criterion

% choosing a function

% PF = @PAL_Weibull; % lin
PF = @PAL_Quick; % lin
% PF = @PAL_CumulativeNormal;
% PF = @PAL_Gumbel;
% PF = @PAL_HyperbolicSecant;

% defining threshold and slope as free parameters and guessing rate and
% lapses as fixed values
paramsFree = [1 1 0 0];

% structure with vector fields .alpha, .beta, .gamma, .lambda collectively 
% defining a 4D parameter grid
% These initial values will serve as seeds for the Nelder-Mead iterative search. 
searchGrid.alpha = [1:0.5:60]; % [0.01:0.001:0.11] [10:.5:60] [1:.05:51]; 
searchGrid.beta = logspace(0,3,119); % logspace(0,3,101)
searchGrid.gamma = .5;
searchGrid.lambda = .02;% 

% uses Nelder-Mead Simplex method to find the maximum in the likelihood function. 
[paramsValues LL exitflag] = PAL_PFML_Fit(StimLevels, NumPos,...
    OutOfNum,searchGrid, paramsFree,PF) % 'checkLimits' 

thres(subi, condi) = paramsValues(1);
exitFlags(subi, condi) = exitflag; 

PropCorrectData = NumPos./OutOfNum;
StimLevelsFine = [min(StimLevels):(max(StimLevels)- ...
min(StimLevels))./1000:max(StimLevels)];
Fit = PF(paramsValues, StimLevelsFine);
plot(StimLevelsFine,Fit,cols{condi},'linewidth',2);
% hold on;
if plotting > 1
%     plot(StimLevels, PropCorrectData, colsp{condi},'markersize',40);
end
set(gca, 'fontsize',12);
% axis([.01 .5 .4 1]);
axis([0 60 .4 1]);


%% Estimating the standard errors
%B = 400;


%[SD paramsSim LLSim converged] = PAL_PFML_BootstrapParametric(StimLevels, OutOfNum, paramsValues, ...
%paramsFree, B, PF, 'searchGrid', searchGrid, 'checkLimits', 0);

% SDs(condi, :) = SD;
% 3.6114    0.2094         0    0.0000
% 1.8264    0.1624         0    0.0000
% 3.8232    0.2485         0    0.0000
% 1.0531    0.1947         0    0.0000

%%

% paramsFree = [1 1 0 1];
% [SD paramsSim LLSim converged] = PAL_PFML_BootstrapNonParametric(StimLevels, NumPos, OutOfNum, ...
% [],paramsFree, B, PF, 'lapseLimits',[0 .03], 'searchGrid', searchGrid);


%% goodness of fit

B = 1000;
[Dev pDev DevSim converged] = PAL_PFML_GoodnessOfFit(StimLevels,...
NumPos, OutOfNum, paramsValues, paramsFree, B, PF, 'searchGrid',...
searchGrid,'checkLimits', 0);

pDevs(condi) = pDev; % 0.4430    0.4430    0.5420    0.2930
devSims(condi,:) = DevSim;
% pause

%% 
end

ylabel('Proportion of correct'); xlabel('Contrast')  
% legend({'Cued High'; ''; 'Cued Low'; ''; 'Non-Cued (high)';'';  'Non-Cued (low)'; ''}) %   24.9548   14.4735   27.6536   15.9465 (26.7926   13.7415   31.4971   16.3844)
if plotting < 2
legend({'Cued High'; 'Cued Low'; 'Non-Cued (high)'; 'Non-Cued (low)'}) %   24.9548   14.4735   27.6536   15.9465 (26.7926   13.7415   31.4971   16.3844)
end
% saveas(figure(subi),[pictureBase '\sub' num2str(subi)], 'jpg');

% end
hold off
%end
% 22.8946   12.5468   27.3640   15.4593
% 0.3170    0.0290    0.7680    0.3490
%% PALAMEDES tutorial

addpath('C:\Users\Richard Naar\Documents\MATLAB\palamedes1_10_3\Palamedes')
% y = PAL_[NameOfFunction] (paramValues, x);

StimLevels = [1:1:6];
pcorrect = PAL_Logistic([3 1 0.5 0],StimLevels); % alfa beta gamma (1/x(-FCT)) lambda 

plot(StimLevels, pcorrect, 'ko');

pcorrect = PAL_Logistic([3 1 0.5],StimLevels);
StimLevels = PAL_Logistic([0 1 0.5 0.01],[0.7 0.8 0.9],'inverse');
FirstDerivs = PAL_Logistic([0 1 0.5 0.01],[-1 0 1 2],'derivative');

% Maximum likelihood criterion

StimLevels = [.01 .03 .05 .07 .09 .11]; 
NumPos = [59 53 68 83 92 99]; 
OutOfNum = [100 100 100 100 100 100];
PF = @PAL_Logistic;

paramsFree = [1 1 0 0];

searchGrid.alpha = [0.01:0.001:0.11];
searchGrid.beta = logspace(0,3,101);
searchGrid.gamma = 0.5;
searchGrid.lambda = 0; % 0.02;

[paramsValues LL exitflag] = PAL_PFML_Fit(StimLevels, NumPos,...
    OutOfNum,searchGrid, paramsFree,PF)

PropCorrectData = NumPos./OutOfNum;
StimLevelsFine = [min(StimLevels):(max(StimLevels)- ...
min(StimLevels))./1000:max(StimLevels)];
Fit = PF(paramsValues, StimLevelsFine);
plot(StimLevelsFine,Fit,'g-','linewidth',2);
hold on;
plot(StimLevels, PropCorrectData,'k.','markersize',40);
set(gca, 'fontsize',12);
axis([0 .12 .4 1]);

% Estimating the standard errors
B = 400;

[SD paramsSim LLSim converged] = PAL_PFML_BootstrapParametric (StimLevels, OutOfNum, paramsValues, ...
paramsFree, B, PF, 'searchGrid', searchGrid);

%% goodness of fit




%% Other
%% load behavioural data
behDir = 'C:\Users\Richard Naar\Documents\dok\ssvep\Visit to York\EEG data\behavioural\';  % directory for datafiles
addpath(behDir)
files = dir([behDir , '*.csv']);

filename = fullfile(behDir,files(1).name);

behDat = readtable(filename,...
    'Delimiter',',','ReadVariableNames',false);

names = table2array(behDat(1,4:41));
names = strrep(names, '.',''); names = strrep(names, '_',''); names = strrep(names, ' ','');

behDat = behDat(3:end,4:41);
behDat.Properties.VariableNames = names;

condi = find(ismember(behDat.trials2label, 'low'));
behDat.trials2response(condi(1))

mean(str2num(cell2mat(behDat.trials2response)))

% behArray = table2array(behDat(:,3:41));
% 
% names = behArray(1,:);
% behArray1 = behArray(3:end,:);
%%