<<<<<<< HEAD
%%  EEG at York prelim (28.11.19)
clc; close all; clear all; 

%% 
do = struct('other',              0, 'plotting',     0,...
            'gedvecs',            0, 'saveCsv',      0,...
            'saveCsvSingle',      0, 'cleanData',    0,...
            'ress',               1);

%% Paths

addpath('C:\Users\Richard Naar\Documents\Matlab toolboxes\ANTeepimport1.13\'); % data folder
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

while do.cleanData == 1
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

%% filter 
lowcutoff = 1;
revfilt = 0; % [0|1] reverse filter (i.e. bandpass filter to notch filter). {default 0}
plotfreqz = 0;
minphase = false; % defalut

[EEG, com, b] = pop_eegfiltnew(EEG, 1, 40, [], revfilt, [], plotfreqz, minphase);
% eegplot(EEG.data, 'events', EEG.event) % see the raw file

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
% juku = eeglab2fieldtrip(EEG,'preprocessing');
%% save the data

save([dataDirClean, implist(subi).name(1:end-4), '-clean.mat'], 'EEG', '-v7.3');
fprintf(['Saveing the dataset as: ' implist(subi).name(1:end-4), '-clean.mat \n'])


end

end

% eegplot(EEG.data, 'events', EEG.event)
%% epoch and save single trial and grand average SNR data

%grandAverage = [];
%p = struct();
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
%% extractin the blink structure
% 
blinkList = dir([dataDir 'blinks\', '*.mat']); 
% load([dataDir 'blinks\' blinkList(end).name])
% fprintf('loading dataset: %s', blinkList(end).name); fprintf('\n') 

for bFileCounter = 1:length(blinkList)
    if strcmp( blinkList( bFileCounter ).name(1:7) , implistM(subi).name(1:7) ) ~= 1 && bFileCounter <= length(blinkList)
        bFileCounter = bFileCounter + 1;
    elseif strcmp( blinkList(bFileCounter).name(1:7) , implistM(subi).name(1:7) ) == 1
        load([dataDir 'blinks\' blinkList(bFileCounter).name])
        fprintf('loading dataset: %s', blinkList(bFileCounter).name); fprintf('\n') % 
    else
        fprintf('Blinker file not found')
    end
end

%% epoch parameters

stimLocked = 1; % stim locked events allow to preserve data near the stim event

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
predDurs = ( [ 2, 2.5, 3, 3.5, 4, 5, 6, 7, 8 ]-transient ) * srate;
% find maximum duration in samples that will be divisible by sampling rate
% (this is needed for coherence averaging)
predDurs = predDurs - mod(predDurs, srate);

% number of trials to skip in each event category
fromTrial = 0; 

%% save same data for the p structure

p.subject(subi).conditionLabels  = event;
% p.conditionLabels = event;
p.transient = transient;
p.fromTrial = fromTrial;
p.subject(subi).psychoLabels = behDatConds; 
p.subject(subi).Id = implistM(subi).name;

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
                eventOnset(currentIndex) = EEG.event(thisEventIndx-1).latency; % find previous event oncet (start if the epoch)
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

electrodes = {'O1', 'Oz' 'O2', 'POz', 'P1', 'P2', 'PO5', 'PO3', 'PO4', 'PO6', 'HEOG'}; % HEOG will be removed latere (NB! hard coded)

p.electrodes = electrodes(1:end-1); % save to p

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
trialId = zeros(65,1);% []; % valid coherence trials
validRess = zeros(65,1); % valid ress trials
rtPsychoTrials = []; rtPsychoTrials = []; correctPsychoTrials = []; % not rally necessary but ok 
%%
while (thisCueEvent <= currentIndex-2) % for thisCueEvent = 1:(currentIndex-1)
     thisCueEvent = thisCueEvent + 1;
     
     startSample = round(eventOnset(thisCueEvent)+transient*srate); % *EEG.srate/1000 % Compensate for the fact that latencies are in ms but data are in samples. Is that really the case?   
     % +(EEG.srate/2)
     endSample = round(eventOffset(thisCueEvent)); % 
     
     trialDur = predDurs(dsearchn(predDurs', endSample-startSample)); % find the trial duration matching the current trial

     % redifine the end or the start of the epoch to be suitable for
     % coherent averaging 
     if stimLocked == 1
         startSample = endSample-trialDur; % find the event oncet if the epoch is calculated from stim event       
     else
         endSample = startSample+trialDur; % find the event offset if the epoch is calculated from cue event
     end 
          
    % Pull out a timeseries that follows/preceeds the event 

    allSamples = EEG.data(elec2plot,startSample:endSample-1); % channels
    
%% prepare ress data matrix
    
    allSamplesRess = EEG.data(1:end-1,startSample:endSample-1); % channels % not including the ocular channel (end-1) 
    rebinnedDataRess = reshape(allSamplesRess, size(allSamplesRess,1), EEG.srate,trialDur/EEG.srate);

%     ressThisEvent = rebinnedDataRe; % not including the ocular channel (end-1) 

%     [artefacts, ~] = findManual(rebinnedDataRess, manualCheck, size(rebinnedDataRess,3), escapeKey, 1);
    artefacts = ones(size(rebinnedDataRess,3),1); % accept all for RESS

    if sum(artefacts) ~= 0 
        validRess(thisCueEvent) = 1;
        ressData{thisCueEvent} = rebinnedDataRess(:,:,find(artefacts)); % there is something worng with it
        ressDataAvg{thisCueEvent} = mean(rebinnedDataRess(:,:,find(artefacts)),3);
    else
    end 
    
%% just to inspect the single trials    
rt = str2num(behDat.trueRT{condi(thisCueEvent+fromTrial)});
isCorrect = behDat.trials2response(condi(thisCueEvent+fromTrial));
contrast = str2num(behDat.trials2intensity{condi(thisCueEvent+fromTrial)});

subtrMean = bsxfun(@minus, allSamples, mean(allSamples,2)); % 
zscored = bsxfun(@rdivide, subtrMean, std(allSamples,0,2)); % 

%%
figure(1)
subplot(2,1,1)
plot(zscored')
hold on

%% find and draw blinks for this trial
 
startSampleInSeconds = startSample/srate;
endSampleInSeconds = endSample/srate;


% extend the window size by this amount 
extS = 0.01;     
extE = 0.01;
if (startSampleInSeconds  - extS) > 0
    startSampleInSeconds = startSampleInSeconds - extS;
end
    endSampleInSeconds = endSampleInSeconds + extE;

% for debugiging 
% % this is just for now start
% endSample = endSampleInSeconds*srate;
% startSample = startSampleInSeconds*srate;
% % jus for now end
    
    
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
correctPsychoTrials(thisCueEvent+fromTrial) = str2num(isCorrect{:});


 
if isempty(rt)
    rtPsychoTrials(thisCueEvent+fromTrial) = nan;
else
    rtPsychoTrials(thisCueEvent+fromTrial) = rt;
end

contrastPsychoTrials(thisCueEvent+fromTrial) = contrast;

%% find a vector containing boolean values for segment rejection (1 keep the trial)
if ii == size(bIndx,2)
    
a = ceil( blinkPeaks/srate);
b = ceil( (blinkPeaks+dur) /srate);
c = ceil( (blinkPeaks-dur) /srate) ;

a(a == 0) = []; b(b == 0) = []; c(c == 0) = [];
a(a > nrSegments) = []; b(b > nrSegments) = []; c(c > nrSegments) = [];

seg = zeros(nrSegments,3);
     
% if size(bIndx,2) == 1
% seg = a+b+c; 
% else    
seg(a,1) = 1;
seg(b,2) = 1;
seg(c,3) = 1;
     
seg = sum(seg,2) ;
% end

notblinks = seg < 2;

end
% txt = find(notblinks == 0)

end

hold off

%% prepare fft data   
    if size(elec2plot,2) > 1
        allSamples=squeeze(mean(allSamples(1:end-1, :))); % mean of channels (not including ocular channel)
    end
    
    % We can now bin into 1 second intervals
    rebinnedData=reshape(allSamples, EEG.srate,trialDur/EEG.srate);
    fftRebinned=fft(rebinnedData); % Perform FFT down time
    % Now the magic happens :) We can average coherently across bins to remove
    % noise
%%
goBack = 0;
[~, goBack] = findManual(rebinnedData, manualCheck, trialDur/EEG.srate, escapeKey, 2);
% notblinks is a vector where 0 indicates a blink

if sum(notblinks) ~= 0 % If there is at least one element that is not a blink
    rebinnedData = rebinnedData(:,find(notblinks));
%     blinkIndex = blinkIndex + 1;

    meanComplexFFT(:,thisCueEvent) = mean(fftRebinned,2); % nb abs % This is the mean complex FFT for this trial (averaged across bins)
%     trialId(thisCueEvent) = thisCueEvent+fromTrial;
        trialId(thisCueEvent+fromTrial) = 1;
    %     meanComplexFFTallTrials(:,blinkIndex) = mean(fftRebinned,2); % keep this for all the trials

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
grandAverage{subi, eventIndx} = mean(meanComplexFFT(:,find(trialId)),2); % Average across trials - the abs means we do NOT still keep coherent information
if ~isempty(meanComplexFFT)   
    allTheSingleTrials{subi, eventIndx} =  abs(meanComplexFFT(1:60,find(trialId)) ); % Here we lose coherence.
    p.allTheSingleTrialsStructure = {'participant', 'condition', 'frequency', 'trials'};
end
allBlinks{subi, eventIndx} = thisEventBlinks;

% for the p
p.subject(subi).condition(eventIndx).data.EEG.artefacts = allBlinks; 
p.subject(subi).condition(eventIndx).data.EEG.validTrials = {find(trialId)}; 
p.subject(subi).condition(eventIndx).data.EEG.single.amplitude = allTheSingleTrials; 
p.subject(subi).condition(eventIndx).data.EEG.average.amplitude = grandAverage; 

p.subject(subi).condition(eventIndx).data.EEG.ressDataSegments = ressData;% squeeze( ressData(subi, eventIndx, find(validRess)) ); 
p.subject(subi).condition(eventIndx).data.EEG.ressDataTrials = ressDataAvg;
p.ressDataStructure = {'participant', 'condition', 'trials', 'time'};

ressAllSegments.subject(subi).condition(eventIndx).data = ressData;% squeeze( ressData(subi, eventIndx, find(validRess)) ); 
ressAllTrials.subject(subi).condition(eventIndx).data = ressDataAvg;% 

% rt, intesity, correct
p.subject(subi).condition(eventIndx).data.psychophysics.correct = {correctPsychoTrials}; 
p.subject(subi).condition(eventIndx).data.psychophysics.rt = {rtPsychoTrials};
p.subject(subi).condition(eventIndx).data.psychophysics.contrast = {contrastPsychoTrials};

%% SNR

hz = abs(grandAverage{subi, eventIndx}(2:44));
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
% %     subplot(2,2,eventIndx); hold on
% %     title(strrep(event{eventIndx}, '_',' '))

%    bar(abs(grandAverage(2:60)));

%     bar(abs(grandAverage{eventIndx}(2:60)));
allSnrE{subi,eventIndx} = snrE;
p.subject(subi).condition(eventIndx).data.EEG.average.snr = snrE; 

% %     bar(allSnrE{subi,eventIndx});

%     stem(abs(grandAverage{eventIndx}(2:60)),'k','linew',3,'markersize',2.5,'markerfacecolor','r')
%     bar(squeeze( mean(snrE(:, eventIndx, 1:end),1)));
end % next event
p.subject(subi).psychoPyTable = behDat; 
% pid = p.subject(subi);
save([[dataDir 'individual\'] date p.subject(subi).Id '.mat'], 'pid', '-v7.3');

end % next subject

% save data ressData
if stimLocked == 1
    save([dataDir date '-TFUR-SNRs-Occipital-stimLocked.mat'], 'allSnrE', '-v7.3');
    save([dataDir date '-TFUR-amp-Occipital-stimLocked.mat'], 'allTheSingleTrials', '-v7.3');
    save([dataDir date '-TFUR-grandAverage-Occipital-stimLocked.mat'], 'grandAverage', '-v7.3');
else
    save([dataDir date '-TFUR-SNRs-Occipital-perdLocked.mat'], 'allSnrE', '-v7.3');
    save([dataDir date '-TFUR-amp-Occipital-perdLocked.mat'], 'allTheSingleTrials', '-v7.3');
end

save([dataDir date '-TFUR-ressAllSegments.mat'], 'ressAllSegments', '-v7.3');
save([dataDir date '-TFUR-ressAllTrials.mat'], 'ressAllTrials', '-v7.3');
save([dataDir date '-allData.mat'], 'p', '-v7.3');

%% blinkStructure
% indiDir = 'C:\Users\Richard Naar\Documents\dok\ssvep\Visit to York\EEG data\individual\'; % data after epoching and artefact inspection
dataList = dir([dataDir, '*.mat']); 

load([dataDir indiList(13).name])
fprintf('Loading: %s \n', indiList(13).name)

bs = struct();

for subi = 1:size(p.subject,2)
    bs.subject(subi).Id = p.subject(subi).Id;
    for condi = 1:length(p.subject(subi).condition)
        bs.subject(subi).condi(1).blinks = p.subject(subi).condition(4).data.EEG.artefacts;
    end
end

save([dataDir date '-blinkStructure.mat'], 'bs', '-v7.3');

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

dNr = 12;
matList = dir([dataDir, '*.mat']); 
load([dataDir matList(12).name])

fprintf('loading dataset: %s', matList(dNr).name); fprintf('\n') % loading dataset: 24-Nov-2019-TFUR-ress.mat

ressAll = ressAllTrials;% ressAllSegments; 
clear snrRESS
for subi = 1:size(ressAll.subject,2) % subjects loop
for ressi = 1:size(ressAll.subject(subi).condition,2) % condition loop
    
data_raw = ressAll.subject(subi).condition(ressi).data';
data_raw = data_raw(~cellfun('isempty',data_raw)); % in case some trials got rejected

cellDat = {}; data = []; dati = 1;
for ij = 1:size(data_raw,1)
    for uj = 1:size(data_raw{ij},3)
        cellDat(dati,1) = {data_raw{ij}(:,:,uj)};
%        data(:,:, dati) = cellDat{dati,1};
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
meanRESS1 = 0; meanRESS2 = 0; meanRESS3 = 0; meanRESS4 = 0; meanRESS5 = 0;
for subi = 1:23
    freso = 0.1;
    f2show = 40;
    f = 2;
    meanRESS1 = meanRESS1 + snrRESS{subi,1,4}(1:(f2show/freso));
    meanRESS2 = meanRESS2 + snrRESS{subi,2,4}(1:(f2show/freso)); 
    meanRESS3 = meanRESS3 + snrRESS{subi,3,4}(1:(f2show/freso)); 
    meanRESS4 = meanRESS4 + snrRESS{subi,4,4}(1:(f2show/freso)); 
%     meanRESS5 = meanRESS5 + snrRESS{subi,1,2}(1:(f2show/freso)); 
    
end
hz2plot = hzR(1:(f2show/freso));

stem(hz2plot, meanRESS1/23, 'k')
hold
stem(hz2plot-.2, meanRESS2/23, 'r')
stem(hz2plot+.2, meanRESS3/23, 'g')
stem(hz2plot+.4, meanRESS4/23, 'b')
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
dati = 1;
load([dataDir matList(dati).name])

fprintf('loading dataset: %s', matList(dati).name); fprintf('\n')

%allSnrE = snrRESS;

iLow = [4:4:50];
iHigh = [15:15:50];
imf = [11, 19, 22, 27, 38, 41];
imfTxt = [{'F2 - F1'},{'F1 + F2'},{'(2*F2) - (2*F1)'},{'3*F1 + F2'},{'(F1 + F2)*2'}, {'3*F2 - F1'}];
allFrex = [iLow, iHigh, imf];

% FOIsLow = zeros(size(subtraction)); FOIsLow(iLow) = subtraction(iLow);

header = [{'SubId', 'cond',...
    'low4', 'low8', 'low12', 'low16', 'low20', 'low24','low28', 'low32', 'low36', 'low40', 'low44', 'low48',...
    'high15', 'high30', 'high45',...
    'F11', 'F19', 'F22','F27', 'F38','F41'}];
data = {};
data(1, :) = header;
subDat = [];
for rind = 1:size(allSnrE,1)
  high = allSnrE{rind,1}(allFrex); % PRED_HIGH
  low =  allSnrE{rind,2}(allFrex); % PRED_LOW
  high_rnd = allSnrE{rind,3}(allFrex); % PRED_RND_H
  low_rnd = allSnrE{rind,4}(allFrex); % PRED_RND_L
  id = num2str(rind);
  
  if rind < 2
      addVal1 = 1; addVal2 = 4;
  else
      addVal1 = addVal1+3 ; addVal2 = addVal2+3;
  end
  
  data(rind+addVal1:rind+addVal2, :) =  [{id}, {'Cued high'}, num2cell(high); {id}, {'Cued low'}, num2cell(low); {id},...
      {'Non-cued high'}, num2cell(high_rnd); {id}, {'Non-cued low'}, num2cell(low_rnd)];
   
end

addpath('C:\Users\Richard Naar\Documents\Matlab toolboxes\letswave6-master\external')

cd(dataDir)
cell2csv([date '-TFUR-SNRs-occipital-stim.csv'], data)

do.saveCsv = 0;
end

%% put the data into table 2 (snrRESS)

while do.saveCsvRess == 1
matList = dir([dataDir, '*.mat']); 
dati = 2;
load([dataDir matList(dati).name])

fprintf('loading dataset: %s', matList(dati).name); fprintf('\n')

iLow = [4:4:50];
iHigh = [15:15:50];
imf = [11, 19, 22, 27, 38, 41];
imfTxt = [{'F2 - F1'},{'F1 + F2'},{'(2*F2) - (2*F1)'},{'3*F1 + F2'},{'(F1 + F2)*2'}, {'3*F2 - F1'}];
allFrex = [iLow, iHigh, imf]/0.1;
snrDat = snrRESS;
% FOIsLow = zeros(size(subtraction)); FOIsLow(iLow) = subtraction(iLow);

header = [{'SubId', 'cond',...
    'low4', 'low8', 'low12', 'low16', 'low20', 'low24','low28', 'low32', 'low36', 'low40', 'low44', 'low48',...
    'high15', 'high30', 'high45',...
    'F11', 'F19', 'F22','F27', 'F38','F41'}];
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

if rind < 2
    addVal1 = 1; addVal2 = 4;
else
    addVal1 = addVal1+3 ; addVal2 = addVal2+3;
end
  

data(rind+addVal1:rind+addVal2, :) =  [{id}, {'Cued high'}, num2cell(high); {id}, {'Cued low'}, num2cell(low); {id},...
      {'Non-cued high'}, num2cell(high_rnd); {id}, {'Non-cued low'}, num2cell(low_rnd)];


end

addpath('C:\Users\Richard Naar\Documents\Matlab toolboxes\letswave6-master\external')

cd(dataDir)
cell2csv([date '-TFUR-SNRs-occipital-stim.csv'], data)

do.saveCsvRess = 0;
end

%% put the data into table

while do.saveCsvSingle == 1
matList = dir([dataDir, '*.mat']); 
load([dataDir matList(2).name])

fprintf('loading dataset: %s', matList(1).name); fprintf('\n')

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
    
  nTrials = size(allTheSingleTrials{rind,1},2);
    
  high = allTheSingleTrials{rind,1}(allFrex,:); % PRED_HIGH
  low =  allTheSingleTrials{rind,2}(allFrex,:); % PRED_LOW
  high_rnd = allTheSingleTrials{rind,3}(allFrex,:); % PRED_RND_H
  low_rnd = allTheSingleTrials{rind,4}(allFrex,:); % PRED_RND_L
  
% for trial = 1:size(low_rnd,2)
  
  id = num2str(rind);
  
  if rind < 2
      addVal1 = 1; addVal2 = 4*size(low_rnd,2);% 4;
  else
      addVal1 = (addVal1+size(low_rnd,2)*4)-1 ; addVal2 = (addVal2+size(low_rnd,2)*4)-1;
  end
  
  data(rind+addVal1:rind+addVal2, :) =  [repmat([{1}, {'Cued high'}],nTrials ,1), num2cell(high)', num2cell(1:nTrials)'; repmat([{id}, {'Cued low'}],nTrials ,1), num2cell(low)', num2cell(1:nTrials)'; ...
      repmat([{id}, {'Non-cued high'}],nTrials ,1), num2cell(high_rnd)', num2cell(1:nTrials)'; repmat([{id}, {'Non-cued low'}],nTrials ,1), num2cell(low_rnd)', num2cell(1:nTrials)'];
  
% end
end

addpath('C:\Users\Richard Naar\Documents\Matlab toolboxes\letswave6-master\external')

cd(dataDir)
cell2csv([date '-TFUR-AllTheSingleTrials-apl-STIM.csv'], data)

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

load([dataDir, implistSNR(8).name]); % load data

fprintf('loading participant: %s \n', implistSNR(8).name);

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
bar(mean( squeeze(cat(3,allSnrE{:,comp1})), 2));
imf = [11, 19, 22, 27, 38, 41];
% imfTxt = ['F2-F1','F1+F2','2*F2-2*F1','3*F1+F2','(F1+F2)*2'];
imfTxt = [{'F2 - F1'},{'F1 + F2'},{'(2*F2) - (2*F1)'},{'3*F1 + F2'},{'(F1 + F2)*2'}, {'3*F2 - F1'}];
maxY = max(mean( squeeze(cat(3,allSnrE{:,comp1})), 2))+1;
% text((imf)-1, zeros(1,length(imf(1:length(imf))))+maxY, imfTxt(1:length(imf))) % *mean(EEG.data(31,:))
set(gca,'ylim',[min(mean( squeeze(cat(3,allSnrE{:,comp1})), 2)) maxY+2], 'FontSize',12)
%% plot the colors
% FOIsLow = zeros(size(grandAverage{1})); FOIsLow(foiLow) = grandAverage{1}(foiLow);
% bar(FOIsLow(2:60), 'r')
% % color high
% FOIsHigh = zeros(size(grandAverage{1})); FOIsHigh(foiHigh) = grandAverage{1}(foiHigh);
%% bar(FOIsHigh(2:60), 'g')
title(strrep(event{comp1}, '_',' '))
%% high - random
% subtraction = grandAverage{comp1}(1:60)-(grandAverage{3}(1:60)+grandAverage{4}(1:60))/2;
subtraction_lrnd = mean( squeeze(cat(3,allSnrE{:,comp1})), 2)-(mean( squeeze(cat(3,allSnrE{:,3})), 2)+mean( squeeze(cat(3,allSnrE{:,4})), 2))/2;

subplot(3,1,2); hold on
bar(subtraction_lrnd);
maxY = max(subtraction_lrnd)+1;
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
set(gca,'ylim',[-4 4], 'FontSize',12)

%% random high
subplot(3,1,3);
avgRnd = (mean( squeeze(cat(3,allSnrE{:,3})), 2)+mean( squeeze(cat(3,allSnrE{:,4})), 2))/2;
bar(avgRnd);
maxY = max(avgRnd)+1;
% text((imf)-1, zeros(1,length(imf(1:length(imf))))+maxY, imfTxt(1:length(imf)))
set(gca,'ylim',[min(avgRnd) maxY+2], 'FontSize',12)
%  title(strrep(event{comp2}, '_',' '))
title('RANDOM')
pause()
end
%% both on the same
% collatz = [30, 15, 46, 23, 70, 35, 106, 53, 160, 80, 40, 20, 10, 5, 16, 8, 4, 2, 1]
comp1 = 1; 
maxF = 40;

% color
colHigh = [0 0 .7];
colLow = [.7 0 0];
imfCol = [0 .7 0];
colPins = [0 0 0];

figN = 1;
figure(figN)
foiLow = [4:4:maxF];  
foiHigh = [15:15:maxF];

% higher frequency
figure(1)
subplot(2,1,1); hold on
hold on
highDat = squeeze(cat(3,allSnrE{:,comp1})); highDat = highDat(1:maxF,:);
randoDat = (squeeze(cat(3,allSnrE{:,3}))+squeeze(cat(3,allSnrE{:,4}))) / 2; randoDat = randoDat(1:maxF,:);
% randoDat = squeeze(cat(3,allSnrE{:,comp1+1})); randoDat = randoDat(1:maxF,:);

subtraction_hrnd = mean( highDat, 2) - mean(randoDat, 2);

bar(subtraction_hrnd, 'FaceColor', colPins);
imf = [11, 19, 22, 27, 38, 41];
imf = imf(find(imf < maxF));

% imfTxt = [{'F2 - F1'},{'F1 + F2'},{'(2*F2) - (2*F1)'},{'3*F1 + F2'},{'(F1 + F2)*2'}, {'3*F2 - F1'}];
% maxY = max(subtraction_hl)+1;
% text((imf)-1, zeros(1,length(imf(1:length(imf))))+maxY, imfTxt(1:length(imf))) % *mean(EEG.data(31,:))

%% plot the colors
% color low

FOIsLow = zeros(size(subtraction_hrnd)); FOIsLow(foiLow) = subtraction_hrnd(foiLow);
bar(FOIsLow, 'FaceColor', colLow)

% color high
FOIsHigh = zeros(size(subtraction_hrnd)); FOIsHigh(foiHigh) = subtraction_hrnd(foiHigh);
bar(FOIsHigh, 'FaceColor', colHigh)
% title([strrep(event{comp1+2}, '_',' '), ' - ',strrep(event{comp1+3}, '_',' ')])

% color imf 
FOIsImf = zeros(size(subtraction_hrnd)); FOIsImf(imf) = subtraction_hrnd(imf);
bar(FOIsImf, 'FaceColor', imfCol)

legend({'Other';'Low';'High';'Imf'})

% bar(FOIsHigh(2:60), 'g')
% title([strrep(event{comp1}, '_',' '), ' - ',strrep(event{comp1+1}, '_',' ')])
title('Attending high')

% 
% err
% mDat = mean( highDat, 2) - mean(randoDat, 2); 
% 
% scaledStd = ( std( highDat - randoDat  ,0,2) ./ sqrt(22) )  * 1.96;
% errHigh = mDat + scaledStd;
% errLow = mDat - scaledStd;
% %
% k = 0.2; ls = 1.5;
% for pli = 1:length(mDat)
%     plot([pli-k pli+k] ,[errHigh(pli) errHigh(pli)], '-', 'LineWidth',ls,'Color', 'k') % upper vertical
%     plot([pli pli] ,[errHigh(pli) errLow(pli)], '-','LineWidth',ls, 'Color', 'k') % horizontal
%     plot([pli-k pli+k] ,[errLow(pli) errLow(pli)], '-','LineWidth',ls, 'Color', 'k') % lower vertical
% end
% % 
%  set(gca,'ylim',[floor(min(errLow)) ceil(max(errHigh))+1], 'FontSize',12)
 set(gca,'ylim',[-3 3], 'FontSize',12)
 %set(gca,'ylim',[floor(min(errLow)) 10], 'FontSize',12)

set(gca,'xlim',[1 40], 'FontSize',12)
ylabel('SNR'); xlabel('Frequency')  


%% high_rnd - low_rnd
lowDat = squeeze(cat(3,allSnrE{:,2})); lowDat = lowDat(1:maxF,:);
% 
% subtraction_rnd = mean( highDat, 2) - mean(randoDat, 2);


subtraction_lrnd = mean( lowDat, 2)- mean( randoDat, 2);

subplot(2,1,2); hold on
bar(subtraction_lrnd, 'FaceColor', colPins);
maxY = max(subtraction_lrnd)+1;
% text((imf)-1, zeros(1,length(imf(1:length(imf))))+maxY, imfTxt(1:length(imf)))
% set(gca,'ylim',[min(subtraction)-2 maxY+2], 'FontSize',12)
%% plot the colors
% color low

FOIsLow = zeros(size(subtraction_lrnd)); FOIsLow(foiLow) = subtraction_lrnd(foiLow);

bar(FOIsLow, 'FaceColor', colLow)


% color high
FOIsHigh = zeros(size(subtraction_lrnd)); FOIsHigh(foiHigh) = subtraction_lrnd(foiHigh);
bar(FOIsHigh, 'FaceColor', colHigh)

% imfTxt = [{'F2 - F1'},{'F1 + F2'},{'(2*F2) - (2*F1)'},{'3*F1 + F2'},{'(F1 + F2)*2'}, {'3*F2 - F1'}];
% maxY = max(subtraction)+1;
% text((imf)-1, zeros(1,length(imf(1:length(imf))))+maxY, imfTxt(1:length(imf))) % *mean(EEG.data(31,:))

% color imf 
FOIsImf = zeros(size(subtraction_lrnd)); FOIsImf(imf) = subtraction_lrnd(imf);
bar(FOIsImf, 'FaceColor', imfCol)

% legend({'Other';'Low';'High';'Imf'})

% err
% mDat = mean( lowDat, 2) - mean(randoDat, 2); 
% 
% 
% scaledStd = ( std( lowDat - randoDat  ,0,2) ./ sqrt(22) )  * 1.96;
% errHigh = mDat + scaledStd;
% errLow = mDat - scaledStd;
% 
% k = 0.2; ls = 1.5;
% for pli = 1:length(mDat)
%     plot([pli-k pli+k] ,[errHigh(pli) errHigh(pli)], '-', 'LineWidth',ls,'Color', 'k') % upper vertical
%     plot([pli pli] ,[errHigh(pli) errLow(pli)], '-','LineWidth',ls, 'Color', 'k') % horizontal
%     plot([pli-k pli+k] ,[errLow(pli) errLow(pli)], '-','LineWidth',ls, 'Color', 'k') % lower vertical
% end
% set(gca,'ylim',[floor(min(errLow)) ceil(max(errHigh))+1], 'FontSize',12)


set(gca,'xlim',[1 40], 'FontSize',12)
set(gca,'ylim',[-3 3], 'FontSize',12)

ylabel('SNR'); xlabel('Frequency')  
title('Attending low')
% title([strrep(event{comp1+2}, '_',' '), ' - ',strrep(event{comp1+3}, '_',' ')])
%%
% figure(2)
% subtraction = highDat - lowDat;
% bar(subtraction)
% 

highDat = squeeze(cat(3,allSnrE{:,comp1})); highDat = highDat(1:maxF,:);
randoLow = squeeze(cat(3,allSnrE{:,3})); randoLow = randoLow(1:maxF,:);
randoHigh = squeeze(cat(3,allSnrE{:,4})); randoHigh = randoHigh(1:maxF,:);
x = 1:length(highDat);

figure()
hold
stem(x-0.2, mean(highDat,2), '-o', 'Color', 'b')
stem(x, mean(lowDat,2), '-o', 'Color', 'r')
stem(x+0.2, mean(randoLow,2), '-o', 'Color', 'k')
stem(x+0.4,mean(randoHigh,2), '-o', 'Color', 'k')

legend({'Cued High';'Cued Low'; 'Non-Cued';})


end

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
% juku = bs.subject(1).condi.blinks{1,1};
% juku = cell(juku);
% cond1 = cell2mat(juku);
%% load in the level averages
cd('C:\Users\Richard Naar\Documents\dok\ssvep\Visit to York\EEG data')

psyList = dir([dataDir, '*.txt']);  %orderfields(implistSNR, date) 
%load([dataDir, psyList(14).name]); % load data

meanContrast = importfile(psyList(5).name, 2, 27);
meanContrast.Properties.VariableNames = {'rowi' 'participant' 'contrast' 'trials_2label' 'NumPos' 'N'};
[C, IA, IC] = unique(meanContrast.trials_2label);
% conds = {'Cued High' 'Cued Low' 'Non-Cued (high)' 'Non-Cued (low)'};
conds = {meanContrast.trials_2label{IA(2)}, meanContrast.trials_2label{IA(1)}, ...
    meanContrast.trials_2label{IA(3)}, meanContrast.trials_2label{IA(4)}};
%% now data for each subject

frex = [8,30,11];
cols = {'r' 'g' 'b' 'k'};
colsp = {'.r' '.g' '.b' '.k'};
% 
% for subi = 2:2%23;
% % subi = 1;
figure()
hold on
for condi = 1:4;
% condi = 3;
frexi = 1;

%% extract data
%EEG
%dataSnr(subi, condi) = allSnrE{subi,condi}(frex(frexi)); 
% condi = 1;
% psycoph
condi = 1;
subCondTable = meanContrast(strcmp(meanContrast.trials_2label, conds{condi}) , :); %& (meanContrast.participant == subi )

OutOfNum = subCondTable.N;
StimLevels = subCondTable.contrast/100;
NumPos = subCondTable.NumPos;

%%
%StimLevels = [.01 .03 .05 .07 .09 .11];
%NumPos = [59 53 68 83 92 99];
%OutOfNum = [100 100 100 100 100 100];
PF = @PAL_Logistic;

paramsFree = [1 1 0 0];

searchGrid.alpha = [0.01:0.001:0.11];
searchGrid.beta = logspace(0,3,101);
searchGrid.gamma = 0.5;
searchGrid.lambda = 0;

[paramsValues LL exitflag] = PAL_PFML_Fit(StimLevels, NumPos,...
    OutOfNum,searchGrid, paramsFree,PF)

paramsValues

PropCorrectData = NumPos./OutOfNum;
StimLevelsFine = [min(StimLevels):(max(StimLevels)- ...
min(StimLevels))./1000:max(StimLevels)];
Fit = PF(paramsValues, StimLevelsFine);
plot(StimLevelsFine,Fit,cols{condi},'linewidth',2);
% hold on;
plot(StimLevels, PropCorrectData, colsp{condi},'markersize',40);
set(gca, 'fontsize',12);
axis([.05 .6 .4 1]);
pause
%% 
end
legend({'Cued High';'Cued Low'; 'Non-Cued (high)'; 'Non-Cued (low)'})


hold off
% end

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
searchGrid.lambda = 0.02;

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

% goodness of fit




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
%% automatic cleaning algorithm
%EEG = clean_rawdata(EEG, 5, [0.25 0.75], 'off','off','off','off'); % , 0.85, 4, 5, 0.25
%
=======
%%  EEG at York prelim (11.11.19)
clc; close all; clear all; 

%% 
do = struct('other',              0, 'plotting',     0,...
            'gedvecs',            0, 'saveCsv',      0,...
            'saveCsvSingle',      0, 'cleanData',    0,...
            'ress',               0);

%% Paths

addpath('C:\Users\Richard Naar\Documents\Matlab toolboxes\ANTeepimport1.13\'); % data folder
dataDir = 'C:\Users\Richard Naar\Documents\dok\ssvep\Visit to York\EEG data\';  % directory for datafiles
dataDirClean = 'C:\Users\Richard Naar\Documents\dok\ssvep\Visit to York\EEG data\clean\';  % directory for datafiles

impdir =    [dataDir, 'raw\'];     
eeglabDir = ' C:\Program Files\MATLAB\R2014a\toolbox\eeglab13_6_5b\';                    % eeglab program directory
addpath(eeglabDir)
addpath(impdir)

implist = dir([impdir, '*.cnt']);          % make a fresh list of the contents of the raw files directory (input of this loop)
implistM = dir([dataDirClean, '*.mat']);   

%% clean, find triggers, down-sample and save

while do.cleanData == 1
    
for subi = 1:length(implist)
%% load data
EEG = [];
EEG = pop_loadeep_v4([impdir, implist(subi).name]); 
fprintf('loading participant: %s \n', implist(subi).name);
%% channel locations
% EEG.chanlocs = ???
% there is also that option
%  EEG = pop_chanedit(EEG, 'load',{locfile 'filetype' 'autodetect'});     % Edit the channel locations structure of an EEGLAB dataset
%% do the cleaning
EEG = clean_drifts(EEG, [0.25 0.75]); % Removes drifts from the data using a forward-backward high-pass filter.
% Maximum tolerated flatline duration. In seconds. If a channel has a longer flatline than this, it will be considered abnormal. 
EEG = clean_flatlines(EEG, 5);

%% DOWN-SAMPLE 
if EEG.srate == 500
else
    EEG = pop_resample(EEG, 500);
end

% eegplot(EEG.data, 'events', EEG.event) 
EEG_ref = mean(EEG.data(1:end,:),1); % 1:64
EEG.data = bsxfun(@minus, EEG.data(:,:), squeeze(EEG_ref)); % average reference

% OUT_EEG = pop_runica( EEG, 'extended', 4, 'interupt', 'on' ); 

%% find triggers
% import trigger information
[~,~,raw] = xlsread([dataDir, 'ViewPixx triggers - 3.xls']); % 3

%% find the triggers % ver3
trigNum = cell2mat(raw(2:5,2:4));
events = zeros(size(EEG.event));
labels = [];
trialCount = 0;
counter = 0;
% 21 1288
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
    
    if e > 1 && ~isempty(strfind(EEG.event(e).type, EEG.event(e-1).type(1:3)))
        EEG.event(e).type = 'erroneusMarker';
    end
    
 labels{1,e} = label;
end

save([dataDirClean, implist(subi).name(1:end-4), '-clean.mat'], 'EEG', '-v7.3');
end
% eegplot(EEG.data, 'events', EEG.event)
end

% eegplot(EEG.data, 'events', EEG.event)
%% epoch with eeglab
% 
% eventTypeRESS = {{'STIM_HIGH' 'STIM_LOW' 'STIM_RND_H' 'STIM_RND_L'} [-10 2; -10 2;-10 2;-10 2] [-200 0; -200 0;-200 0;-200 0]}; 
% % eventTypeRESS = {{'STIM_HIGH'} [-10 2] [-200 0]}; 
% 
% EEG1 = pop_epoch(EEG, eventTypeRESS{1}, eventTypeRESS{2}, 'epochinfo', 'yes');  % extract epochs 1 second before the picture and 3 after the picture. We may want to shortening
% % EEG = pop_rmbase(EEG, eventTypeRESS{3});                                    % remove baseline
% 
% EEG1 = pop_runica(EEG, 'extended', 1, 'interupt', 'on');     % running the ICA    
%% epoch and save single trial and grand average SNR data

grandAverage = [];

%% loop through all the datasets in implistM (pre-processed data) one by one

for subi = 1:3% length(implistM)
EEG = []; 

load([dataDirClean, implistM(subi).name]); % load data

fprintf('loading participant: %s \n', implistM(subi).name);

srate = EEG.srate;  % sampling frequency 
%% epoch parameters

stimLocked = 1; % stim locked events allow to preserve data near the stim event

if stimLocked == 1
    event = {'STIM_HIGH' 'STIM_LOW' 'STIM_RND_H' 'STIM_RND_L'}; 
else
    event = {'PRED_HIGH' 'PRED_LOW' 'PRED_RND_H' 'PRED_RND_L'};
end

% the epochs will defined based on the endpoint (e.g. 500 samples before
% the stimulus marker 'SRIM_HIGH')
% event = {'STIM_HIGH' 'STIM_LOW' 'STIM_RND_H' 'STIM_RND_L'}; 
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
predDurs = ( [ 2, 2.5, 3, 3.5, 4, 5, 6, 7, 8 ]-transient ) * srate;
% find maximum duration in samples that will be divisible by sampling rate
% (this is needed for coherence averaging)
predDurs = predDurs - mod(predDurs, srate);

% number of trials to skip in each event category
fromTrial = 5; 

%%
for eventIndx = 1:length(event) % loop through the event categories

% Find all the events corresponding to the eventIndex
thisEvent = event{eventIndx}; % pick events from predefined event categories
nEvents=length(EEG.event); % get the number of all the events in the dataset
currentIndex=1; % set the eventindex to be one before entering the loop 
% (this will be increased by one each time we find the predefined event from 
% the eventlist, after exclusion of the number of trials set in freomTrial)
eventCounter = 0; % set the counter to exclude first fwe trials

for thisEventIndx=1:nEvents % loop throgh all the events
    if  strfind(EEG.event(thisEventIndx).type, char(thisEvent)) % any(ismember(thisEvent, EEG.event(thisEventIndx).type)) % to compare cued to non-cued
        eventCounter = eventCounter + 1;
        if eventCounter > fromTrial %
            if stimLocked == 1
                eventOffset(currentIndex) = EEG.event(thisEventIndx).latency; % find event offset
                eventOnset(currentIndex) = EEG.event(thisEventIndx-1).latency; % find previous event oncet (start if the epoch)
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

electrodes = {'O1', 'Oz' 'O2', 'POz', 'P1', 'P2', 'PO5', 'PO3', 'PO4', 'PO6'}; 

elec2plot = find(ismember({EEG.chanlocs.labels}, electrodes));
elec2plotNames = {EEG.chanlocs(elec2plot).labels};

fprintf('\nNumber of electrodes aggregated: %d ', length(elec2plot)); fprintf('\n')
% elec2plotNames 
%% epoch the data
% ressProov = double.empty(4, 64, 3500,0);
for thisCueEvent=1:(currentIndex-1) 

     startSample = round(eventOnset(thisCueEvent)+transient*srate); % *EEG.srate/1000 % Compensate for the fact that latencies are in ms but data are in samples. Is that really the case?   
     % +(EEG.srate/2)
     endSample = round(eventOffset(thisCueEvent)); % 
     
     trialDur = predDurs(dsearchn(predDurs', endSample-startSample)); % find the trial duration matching the current trial
     
     ressData{subi, eventIndx, thisCueEvent, :} = EEG.data(:,startSample:endSample);
     
     % redifine the end or the start of the epoch to be suitable for
     % coherent averaging 
     if stimLocked == 1
         startSample = endSample-trialDur; % find the event oncet if the epoch is calculated from stim event       
     else
         endSample = startSample+trialDur; % find the event offset if the epoch is calculated from cue event
     end 
     
      
          
    % Pull out a timeseries that follows/preceeds the event 

    allSamples=EEG.data(elec2plot,startSample:endSample-1); % channels
    if size(elec2plot,2) > 1
        allSamples=squeeze(mean(allSamples)); % mean of channels
    end
    % We can now bin into 1 second intervals
    rebinnedData=reshape(allSamples,EEG.srate,trialDur/EEG.srate);
    fftRebinned=fft(rebinnedData); % Perform FFT down time
    % Now the magic happens :) We can average coherently across bins to remove
    % noise

    meanComplexFFT(:,thisCueEvent)=mean(fftRebinned,2); % nb abs % This is the mean complex FFT for this trial (averaged across bins)
end

% grandAverage{eventIndx}=mean((meanComplexFFT),2); % Average across trials - the abs means we do NOT still keep coherent information
grandAverage{subi, eventIndx} = mean((meanComplexFFT),2); % Average across trials - the abs means we do NOT still keep coherent information
allTheSingleTrials{subi, eventIndx} = abs(meanComplexFFT(1:50, :));

%% SNR

hz = abs(grandAverage{subi, eventIndx}(2:90));
snrE = zeros(1,size(hz,1));
skipbins =  1; % 1 Hz, hard-coded!
numbins  = 2; %  2 Hz, also hard-coded!

% loop over frequencies and compute SNR
for hzi=numbins+1:length(hz)-numbins-1
    numer = hz(hzi);
    denom = mean( hz([hzi-numbins:hzi-skipbins hzi+skipbins:hzi+numbins]) );
    snrE(hzi) = numer./denom;
end  
%%
    subplot(2,2,eventIndx); hold on
    title(strrep(event{eventIndx}, '_',' '))

%    bar(abs(grandAverage(2:60)));

%     bar(abs(grandAverage{eventIndx}(2:60)));
    allSnrE{subi,eventIndx} = snrE;
    bar(allSnrE{subi,eventIndx});

%     stem(abs(grandAverage{eventIndx}(2:60)),'k','linew',3,'markersize',2.5,'markerfacecolor','r')
%     bar(squeeze( mean(snrE(:, eventIndx, 1:end),1)));
end
end

% save data ressData
if stimLocked == 1
    save([dataDir date '-TFUR-SNRs-Occipital-stimLocked.mat'], 'allSnrE', '-v7.3');
    save([dataDir date '-TFUR-amp-Occipital-stimLocked.mat'], 'allTheSingleTrials', '-v7.3');
else
    save([dataDir date '-TFUR-SNRs-Occipital-perdLocked.mat'], 'allSnrE', '-v7.3');
    save([dataDir date '-TFUR-amp-Occipital-perdLocked.mat'], 'allTheSingleTrials', '-v7.3');
end

save([dataDir date '-TFUR-ress.mat'], 'ressData', '-v7.3');

%% Rhythmic Entrainment Source Separation
% only for lower frequency at the moment
while do.ress
% parameters for RESS:
srate = 500;
peakwidt  = .5; % FWHM at peak frequency
neighfreq = 1;  % distance of neighboring frequencies away from peak frequency, +/- in Hz
neighwidt = 1;  % FWHM of the neighboring frequencies

freq1 = 8;
% freq2 = 30;

nfft = ceil( EEG.srate/.1 ); % .1 Hz resolution
hzR    = linspace(0,EEG.srate,nfft);

dNr = 1;
matList = dir([dataDir, '*.mat']); 
load([dataDir matList(dNr).name])

fprintf('loading dataset: %s', matList(dNr).name); fprintf('\n')

for subi = 1:size(ressData,1)
for ressi = 1:size(ressData,2)    

data = squeeze(ressData(subi,ressi,:))';   

data = cell2mat( data );    

nbchan = size(data,1);
% data = p.data_stacked ;    
% GED for spatial filter

[covHi, covLo, covAt] = deal( zeros(nbchan) );
    
% compute covariance matrix for lower neighbor
fdatLo = filterFGx(data,srate,freq1+neighfreq,neighwidt);
fdatLo = reshape( fdatLo(:,:,:), nbchan,[] );
fdatLo = bsxfun(@minus,fdatLo,mean(fdatLo,2));
covLo  = covLo + (fdatLo*fdatLo')/size(data,2); % pnts?

% compute covariance matrix for upper neighbor
fdatHi = filterFGx(data,srate,freq1-neighfreq,neighwidt);
fdatHi = reshape( fdatHi(:,:,:), nbchan,[] );
fdatHi = bsxfun(@minus,fdatHi,mean(fdatHi,2));
covHi  = covHi + (fdatHi*fdatHi')/size(data,2);
    
% compute covariance matrix at peak frequency
fdatAt = filterFGx(data,srate,freq1,peakwidt);
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

ress_ts{subi, ressi} = ress_ts1;

%% compute SNR spectrum for RESS

ressx = mean(abs( fft(ress_ts{subi, ressi}(:,:),nfft,1)/size(data,2) ).^2,2);

snrR = deal(zeros(size(hzR)));
skipbins =  10; % 1 Hz, hard-coded!
numbins  = 20+skipbins; %  2 Hz, also hard-coded!


% loop over frequencies and compute SNR
for hzi=numbins+1:length(hzR)-numbins-1
    numer = ressx(hzi);
    denom = mean( ressx([hzi-numbins:hzi-skipbins hzi+skipbins:hzi+numbins]) );
    snrR(hzi) = numer./denom;
end

snrRESS{subi, ressi} = snrR;

end
end

do.ress = 0;
end

%%

plot(hzR(1:(40/0.1)), snrRESS{1,1}(1:(40/0.1)), 'k')
hold
plot(hzR(1:(40/0.1)), snrRESS{1,2}(1:(40/0.1)), 'r')
plot(hzR(1:(40/0.1)), snrRESS{1,3}(1:(40/0.1)), 'g')
plot(hzR(1:(40/0.1)), snrRESS{1,4}(1:(40/0.1)), 'b')
hold
%% put the data into table 1

while do.saveCsv == 1
matList = dir([dataDir, '*.mat']); 
load([dataDir matList(1).name])

fprintf('loading dataset: %s', matList(1).name); fprintf('\n')


iLow = [4:4:50];
iHigh = [15:15:50];
imf = [11, 19, 22, 27, 38, 41];
imfTxt = [{'F2 - F1'},{'F1 + F2'},{'(2*F2) - (2*F1)'},{'3*F1 + F2'},{'(F1 + F2)*2'}, {'3*F2 - F1'}];
allFrex = [iLow, iHigh, imf];

addpath('C:\Users\Richard Naar\Documents\fNIRS\Kersteni uurimistöö\data\analysisR\fnirs')

% FOIsLow = zeros(size(subtraction)); FOIsLow(iLow) = subtraction(iLow);

header = [{'SubId', 'cond',...
    'low4', 'low8', 'low12', 'low16', 'low20', 'low24','low28', 'low32', 'low36', 'low40', 'low44', 'low48',...
    'high15', 'high30', 'high45',...
    'F11', 'F19', 'F22','F27', 'F38','F41'}];
data = {};
data(1, :) = header;
subDat = [];
for rind = 1:size(allSnrE,1)
  high = allSnrE{rind,1}(allFrex); % PRED_HIGH
  low =  allSnrE{rind,2}(allFrex); % PRED_LOW
  high_rnd = allSnrE{rind,3}(allFrex); % PRED_RND_H
  low_rnd = allSnrE{rind,4}(allFrex); % PRED_RND_L
  id = num2str(rind);
  
  if rind < 2
      addVal1 = 1; addVal2 = 4;
  else
      addVal1 = addVal1+3 ; addVal2 = addVal2+3;
  end
  
  data(rind+addVal1:rind+addVal2, :) =  [{id}, {'Cued high'}, num2cell(high); {id}, {'Cued low'}, num2cell(low); {id},...
      {'Non-cued high'}, num2cell(high_rnd); {id}, {'Non-cued low'}, num2cell(low_rnd)];
   
end

addpath('C:\Users\Richard Naar\Documents\Matlab toolboxes\letswave6-master\external')

cd(dataDir)
cell2csv([date '-TFUR-SNRs-occipital-stim.csv'], data)

do.saveCsv = 0;
end

%% put the data into table

while do.saveCsvSingle == 1
matList = dir([dataDir, '*.mat']); 
load([dataDir matList(2).name])

fprintf('loading dataset: %s', matList(1).name); fprintf('\n')

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
    
  nTrials = size(allTheSingleTrials{rind,1},2);
    
  high = allTheSingleTrials{rind,1}(allFrex,:); % PRED_HIGH
  low =  allTheSingleTrials{rind,2}(allFrex,:); % PRED_LOW
  high_rnd = allTheSingleTrials{rind,3}(allFrex,:); % PRED_RND_H
  low_rnd = allTheSingleTrials{rind,4}(allFrex,:); % PRED_RND_L
  
% for trial = 1:size(low_rnd,2)
  
  id = num2str(rind);
  
  if rind < 2
      addVal1 = 1; addVal2 = 4*size(low_rnd,2);% 4;
  else
      addVal1 = (addVal1+size(low_rnd,2)*4)-1 ; addVal2 = (addVal2+size(low_rnd,2)*4)-1;
  end
  
  data(rind+addVal1:rind+addVal2, :) =  [repmat([{1}, {'Cued high'}],nTrials ,1), num2cell(high)', num2cell(1:nTrials)'; repmat([{id}, {'Cued low'}],nTrials ,1), num2cell(low)', num2cell(1:nTrials)'; ...
      repmat([{id}, {'Non-cued high'}],nTrials ,1), num2cell(high_rnd)', num2cell(1:nTrials)'; repmat([{id}, {'Non-cued low'}],nTrials ,1), num2cell(low_rnd)', num2cell(1:nTrials)'];
  
% end
end

addpath('C:\Users\Richard Naar\Documents\Matlab toolboxes\letswave6-master\external')

cd(dataDir)
cell2csv([date '-TFUR-AllTheSingleTrials-apl-STIM.csv'], data)

do.saveCsvSingle = 0;
end

% plotting ERP

while do.plotting == 1

    figure(1)
    hold
    col = {'r','g','b','k'};

    for condi = 1:size(grandAverageERP,2)

    conERPd = mean( squeeze(cat(3,grandAverageERP{:,condi})), 2);
    baseline = mean(conERPd(1:500));

    plot(conERPd-baseline, col{condi});
    %plot(smooth(conERPd-baseline, 100), col{condi});

    end
    do.plotting = 0;
    
end
%%
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
[x,idx] = sort(datenum([{implistSNR.date}]), 'descend');
implistSNR = implistSNR(idx);

load([dataDir, implistSNR(1).name]); % load data

fprintf('loading participant: %s \n', implistSNR(1).name);

for frex = 1:2
comp1 = frex; 
% comp2 = 3;
figN = 1;
% figure(figN)
foiLow = [4:4:59];  %foiLow = foiLow+1;
foiHigh = [15:15:59];%  foiHigh = foiHigh+1;
% high
figure(frex)
subplot(3,1,1); hold on
% bar(grandAverage{comp1}(2:60));
bar(mean( squeeze(cat(3,allSnrE{:,comp1})), 2));
imf = [11, 19, 22, 27, 38, 41];
% imfTxt = ['F2-F1','F1+F2','2*F2-2*F1','3*F1+F2','(F1+F2)*2'];
imfTxt = [{'F2 - F1'},{'F1 + F2'},{'(2*F2) - (2*F1)'},{'3*F1 + F2'},{'(F1 + F2)*2'}, {'3*F2 - F1'}];
maxY = max(mean( squeeze(cat(3,allSnrE{:,comp1})), 2))+1;
text((imf)-1, zeros(1,length(imf(1:length(imf))))+maxY, imfTxt(1:length(imf))) % *mean(EEG.data(31,:))
set(gca,'ylim',[min(mean( squeeze(cat(3,allSnrE{:,comp1})), 2)) maxY+2], 'FontSize',12)
%% plot the colors
% FOIsLow = zeros(size(grandAverage{1})); FOIsLow(foiLow) = grandAverage{1}(foiLow);
% bar(FOIsLow(2:60), 'r')
% % color high
% FOIsHigh = zeros(size(grandAverage{1})); FOIsHigh(foiHigh) = grandAverage{1}(foiHigh);
%% bar(FOIsHigh(2:60), 'g')
title(strrep(event{comp1}, '_',' '))
%% high - random
% subtraction = grandAverage{comp1}(1:60)-(grandAverage{3}(1:60)+grandAverage{4}(1:60))/2;
subtraction = mean( squeeze(cat(3,allSnrE{:,comp1})), 2)-(mean( squeeze(cat(3,allSnrE{:,3})), 2)+mean( squeeze(cat(3,allSnrE{:,4})), 2))/2;

subplot(3,1,2); hold on
bar(subtraction);
maxY = max(subtraction)+1;
% text((imf)-1, zeros(1,length(imf(1:length(imf))))+maxY, imfTxt(1:length(imf)))
% set(gca,'ylim',[min(subtraction)-2 maxY+2], 'FontSize',12)
%% plot the colors
% color low

FOIsLow = zeros(size(subtraction)); FOIsLow(foiLow) = subtraction(foiLow);
bar(FOIsLow, 'm')

% color high
FOIsHigh = zeros(size(subtraction)); FOIsHigh(foiHigh) = subtraction(foiHigh);
bar(FOIsHigh, 'k')
title('certain - uncertain condition')
legend({'Other';'Low';'High'})

%% random high
subplot(3,1,3);
avgRnd = (mean( squeeze(cat(3,allSnrE{:,3})), 2)+mean( squeeze(cat(3,allSnrE{:,4})), 2))/2;
bar(avgRnd);
maxY = max(avgRnd)+1;
text((imf)-1, zeros(1,length(imf(1:length(imf))))+maxY, imfTxt(1:length(imf)))
set(gca,'ylim',[min(avgRnd) maxY+2], 'FontSize',12)
%  title(strrep(event{comp2}, '_',' '))
title('RANDOM')
pause()
end
%% both on the same
% collatz = [30, 15, 46, 23, 70, 35, 106, 53, 160, 80, 40, 20, 10, 5, 16, 8, 4, 2, 1]
comp1 = 1; 
% comp2 = 3;
figN = 1;
figure(figN)
foiLow = [4:4:59];  %foiLow = foiLow+1;
foiHigh = [15:15:59];%  foiHigh = foiHigh+1;
% high
figure(1)
subplot(2,1,1); hold on
% bar(grandAverage{comp1}(2:60));
subtraction_hl = mean( squeeze(cat(3,allSnrE{:,comp1})), 2) - mean( squeeze(cat(3,allSnrE{:,comp1+1})), 2);
bar(subtraction_hl);
imf = [11, 19, 22, 27, 38, 41];
% imfTxt = ['F2-F1','F1+F2','2*F2-2*F1','3*F1+F2','(F1+F2)*2'];
imfTxt = [{'F2 - F1'},{'F1 + F2'},{'(2*F2) - (2*F1)'},{'3*F1 + F2'},{'(F1 + F2)*2'}, {'3*F2 - F1'}];
maxY = max(subtraction_hl)+1;
text((imf)-1, zeros(1,length(imf(1:length(imf))))+maxY, imfTxt(1:length(imf))) % *mean(EEG.data(31,:))
set(gca,'ylim',[-3 3], 'FontSize',12)
%% plot the colors
% color low

FOIsLow = zeros(size(subtraction_hl)); FOIsLow(foiLow) = subtraction_hl(foiLow);
bar(FOIsLow, 'm')

% color high
FOIsHigh = zeros(size(subtraction_hl)); FOIsHigh(foiHigh) = subtraction_hl(foiHigh);
bar(FOIsHigh, 'k')
title([strrep(event{comp1+2}, '_',' '), ' - ',strrep(event{comp1+3}, '_',' ')])
legend({'Other';'Low';'High'})

% bar(FOIsHigh(2:60), 'g')
title([strrep(event{comp1}, '_',' '), ' - ',strrep(event{comp1+1}, '_',' ')])
%% high_rnd - low_rnd
% subtraction = grandAverage{comp1}(1:60)-(grandAverage{3}(1:60)+grandAverage{4}(1:60))/2;
subtraction = mean( squeeze(cat(3,allSnrE{:,3})), 2)-mean( squeeze(cat(3,allSnrE{:,4})), 2);

subplot(2,1,2); hold on
bar(subtraction);
maxY = max(subtraction)+1;
% text((imf)-1, zeros(1,length(imf(1:length(imf))))+maxY, imfTxt(1:length(imf)))
% set(gca,'ylim',[min(subtraction)-2 maxY+2], 'FontSize',12)
%% plot the colors
% color low

FOIsLow = zeros(size(subtraction)); FOIsLow(foiLow) = subtraction(foiLow);
bar(FOIsLow, 'm')

% color high
FOIsHigh = zeros(size(subtraction)); FOIsHigh(foiHigh) = subtraction(foiHigh);
bar(FOIsHigh, 'k')

imfTxt = [{'F2 - F1'},{'F1 + F2'},{'(2*F2) - (2*F1)'},{'3*F1 + F2'},{'(F1 + F2)*2'}, {'3*F2 - F1'}];
maxY = max(subtraction)+1;
text((imf)-1, zeros(1,length(imf(1:length(imf))))+maxY, imfTxt(1:length(imf))) % *mean(EEG.data(31,:))

legend({'Other';'Low';'High'})

set(gca,'ylim',[-3 3], 'FontSize',12)
title([strrep(event{comp1+2}, '_',' '), ' - ',strrep(event{comp1+3}, '_',' ')])

end

%% Other
%% load behavioural data
behDir = 'C:\Users\Richard Naar\Documents\dok\ssvep\Visit to York\EEG data\behavioural\';  % directory for datafiles
addpath(behDir)
files = dir([behDir , '*.csv']);

filename = fullfile(behDir,files(2).name);
% strcmp(files(1).name, 'AW') 
behDat = readtable(filename,...
    'Delimiter',',','ReadVariableNames',false);
behArray = table2array(behDat(:,2:41));

%%
>>>>>>> 7714be44f13835ed15402eeb51ca8a0c88470e9d
