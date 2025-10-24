%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data as first time (before saving)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
close all;

region = 'ITC';
cd('E:\zebel_jenab');
face_name = ['Face_Mean_PSTH_data_' , region ,'.mat'];
body_name = ['body_Mean_PSTH_data_' , region ,'.mat'];
natural_name = ['natural_Mean_PSTH_data_' , region ,'.mat'];
artifact_name = ['artifact_Mean_PSTH_data_' , region ,'.mat'];
cd('E:\zebel_jenab');
[SpikeTrain_it_all, folderNames] = load_data_lfp(region);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
number_of_neurons       = length(SpikeTrain_it_all);
window_length           = 50;
sliding_step            = 5;
dt                      = 0.001;
iteration               = 5;

number_of_time_slices   = floor(((900 - window_length)/sliding_step) + 1);

Max_of_confidece_int    = ceil((iteration*97.5)/100);

Min_of_confidece_int    = ceil((iteration*2.5)/100);

min_of_number_trials    = 819;

smooth_width            = 50;
nStimuli                = 144;

trainRatio              = 0.7;
face_labels             = dlmread('face_labels.txt');
body_labels             = dlmread('body_labels.txt');
artifact_labels         = dlmread('artifact_labels.txt');
natural_labels          = dlmread('natural_labels.txt');
nonface_labels          = dlmread('nonface_labels.txt');

body_labels             = body_labels(1:length(artifact_labels));
face_labels             = face_labels(1:length(artifact_labels));

natural_labels          = natural_labels(1:length(artifact_labels));

face_counter            = 0;
face_counter_data       = 0;
nonFCounter             = 0;

min_stimulus            = 50;

for i = 3:length(SpikeTrain_it_all)
    SpikeTrain_it_all(i).data = SpikeTrain_it_all(i).data.ua;
    SpikeTrain_it_all(i).cm = SpikeTrain_it_all(i).cm.cm;
    SpikeTrain_it_all(i).data_it = SpikeTrain_it_all(i).data_it.data;
    SpikeTrain_it_all(i).data_pfc = SpikeTrain_it_all(i).data_pfc.data;
end

time = linspace(-200,700,number_of_time_slices);
face_counter      = 0;
face_counter_data = 0;
nonFCounter       = 0;
mean_vec_face     = zeros(number_of_neurons , number_of_time_slices);
var_vec_face      = zeros(number_of_neurons , number_of_time_slices);
mean_vec_artifact = zeros(number_of_neurons , number_of_time_slices);
var_vec_artifact  = zeros(number_of_neurons , number_of_time_slices);
mean_vec_natural  = zeros(number_of_neurons , number_of_time_slices);
var_vec_natural   = zeros(number_of_neurons , number_of_time_slices);
mean_vec_body     = zeros(number_of_neurons , number_of_time_slices);
var_vec_body      = zeros(number_of_neurons , number_of_time_slices);
mean_vec_nonface  = zeros(number_of_neurons , number_of_time_slices);
var_vec_nonface   = zeros(number_of_neurons , number_of_time_slices);



for i= 3:length(SpikeTrain_it_all)
    
    if (size(SpikeTrain_it_all(i).data_it,2) > 900)
        SpikeTrain_it_all(i) = [];
    end
end

save(['LFP_data','.mat'] , 'SpikeTrain_it_all' ,'-v7.3' );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load data
region = 'ITC';
cd('E:\zebel_jenab');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
SpikeTrain_it_all       = load('LFP_data.mat');
SpikeTrain_it_all       = SpikeTrain_it_all.SpikeTrain_it_all;

number_of_neurons       = length(SpikeTrain_it_all);
window_length           = 50;
sliding_step            = 5;
dt                      = 0.001;
iteration               = 5;

number_of_time_slices   = floor(((900 - window_length)/sliding_step) + 1);

Max_of_confidece_int    = ceil((iteration*97.5)/100);

Min_of_confidece_int    = ceil((iteration*2.5)/100);

min_of_number_trials    = 819;

smooth_width            = 50;
nStimuli                = 144;

trainRatio              = 0.7;
face_labels             = dlmread('face_labels.txt');
body_labels             = dlmread('body_labels.txt');
artifact_labels         = dlmread('artifact_labels.txt');
natural_labels          = dlmread('natural_labels.txt');
nonface_labels          = dlmread('nonface_labels.txt');

body_labels             = body_labels(1:length(artifact_labels));
face_labels             = face_labels(1:length(artifact_labels));

natural_labels          = natural_labels(1:length(artifact_labels));

face_counter            = 0;
face_counter_data       = 0;
nonFCounter             = 0;

min_stimulus            = 50;


time = linspace(-200,700,number_of_time_slices);
face_counter      = 0;
face_counter_data = 0;
nonFCounter       = 0;
mean_vec_face     = zeros(number_of_neurons , number_of_time_slices);
var_vec_face      = zeros(number_of_neurons , number_of_time_slices);
mean_vec_artifact = zeros(number_of_neurons , number_of_time_slices);
var_vec_artifact  = zeros(number_of_neurons , number_of_time_slices);
mean_vec_natural  = zeros(number_of_neurons , number_of_time_slices);
var_vec_natural   = zeros(number_of_neurons , number_of_time_slices);
mean_vec_body     = zeros(number_of_neurons , number_of_time_slices);
var_vec_body      = zeros(number_of_neurons , number_of_time_slices);
mean_vec_nonface  = zeros(number_of_neurons , number_of_time_slices);
var_vec_nonface   = zeros(number_of_neurons , number_of_time_slices);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sep stimuli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
SpikeTrain_it_all = seperateStimulus_plv(SpikeTrain_it_all,face_labels,'face',177);
SpikeTrain_it_all = seperateStimulus_plv(SpikeTrain_it_all,body_labels,'body',177);
SpikeTrain_it_all = seperateStimulus_plv(SpikeTrain_it_all,natural_labels,'natural',177);
SpikeTrain_it_all = seperateStimulus_plv(SpikeTrain_it_all,artifact_labels,'artifact',177);
SpikeTrain_it_all = seperateStimulus_plv(SpikeTrain_it_all,nonface_labels,'nonface',177);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% resize task information (cm) to be equal to number of trials 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i= 3:length(SpikeTrain_it_all)
    SpikeTrain_it_all(i).cm = SpikeTrain_it_all(i).cm(1:end -2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculation neural Varibility 

neuralVar = computeTrialByTrialVariability(SpikeTrain_it_all, window_length, sliding_step);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% cpmputes and rank neural varibilty to 2 time window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;

time_range_early = [60 290]; % early window --> first   quenching on MMFF
time_range_late  = [300 500];% late window  --> Second  quenching on MMFF
topPercent = 0.2;            % Percent of top and down Neural varibility trials that picks for continue process 

[highVarTrials_early, lowVarTrials_early, rankedStim_early] = ...
    findExtremeVariabilityTrials_timeWindow(SpikeTrain_it_all, neuralVar, topPercent, time_range_early, window_length, sliding_step);

[highVarTrials_late, lowVarTrials_late, rankedStim_late] = ...
    findExtremeVariabilityTrials_timeWindow(SpikeTrain_it_all, neuralVar, topPercent, time_range_late, window_length, sliding_step);



[highFR_lowVar_early, lowFR_lowVar_early, rankedRate_early, rankedVar_early] = ...
    findHighAndLowRateLowVarStim(neuralVar, time_range_early, topPercent, window_length, sliding_step);


[highFR_lowVar_late, lowFR_lowVar_late, rankedRate_late, rankedVar_late] = ...
    findHighAndLowRateLowVarStim(neuralVar, time_range_early, topPercent, window_length, sliding_step);


disp('Neuron 1: High-rate & stable stimuli:');
disp(highFR_lowVar_early{3});

disp('Neuron 1: Low-rate & stable stimuli:');
disp(lowFR_lowVar_early{3});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---- Edited PLV section (replace the original PLV loop & plotting) ----
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
freq = [1:1:45 55:1:160];%[4:5:55, 65:10:130];
% Ensure we use the usable sample range the user mentioned


% reuse frequency vector 'freq' already defined earlier
fs = 1000; pad = 1;        % sampling / wavelet settings (as before)

N_shuff = 500;             % number of shuffles for baseline
nConds = 4;                % number of conditions (you used 4)

% Preallocate containers
PLV_it_pfc = nan(nConds, length(freq), 900);
PLV_pfc_it_shuff = nan(nConds, N_shuff, length(freq), 900);
PLV_shuf_c = nan(nConds, length(freq), 900);
phaseDiff = cell(1,nConds);
sig_it_wav = cell(1,nConds);
sig_pfc_wav = cell(1,nConds);

% Extract LFPs and restrict to usable_idx
sig_it = SpikeTrain_it_all(ss).data_it;     % trials x time
sig_pfc = SpikeTrain_it_all(ss).data_pfc;

% Remove artifacts (use your artifact function)
artif_it  = lfp_artifact_remove(SpikeTrain_it_all(ss).data_it, 0);
artif_pfc = lfp_artifact_remove(SpikeTrain_it_all(ss).data_pfc, 0);
artif = artif_it | artif_pfc;

% Compute complex wavelet transform (time window already limited to usable_idx)
[sig_it_an, fr]   = ndass_wavelet_np(sig_it, freq, fs, pad);   % dims: trials x freq x time
[sig_pfc_an, fr]  = ndass_wavelet_np(sig_pfc, freq, fs, pad);  % same dims

L = size(sig_it_an, 3);  % number of time points returned by the wavelet function
% if L ~= length(usable_idx) consider aligning (but usually matches)
%
% Loop conditions
for ci = 1:nConds
    fprintf('Computing PLV for condition %d / %d\n', ci, nConds);
    % find trials for this condition (same logic as your original code)
%     ind_h = ismember(R(ss).TriaInf(:,1), level_occ_4l(:,ci)) & (~artif);
    ind_h = highFR_lowVar_late{ss};
    % get phase matrices: dims -> trials_condition x freq x time
    phi_it = squeeze(angle(sig_it_an(ind_h, :, :)));    % (nTrials_cond x freq x time)
    phi_pfc = squeeze(angle(sig_pfc_an(ind_h, :, :)));

    % save raw wavelets for this condition (if you want)
    sig_it_wav{ci} = sig_it_an(ind_h,:,:);
    sig_pfc_wav{ci} = sig_pfc_an(ind_h,:,:);

    % compute PLV across trials for each freq x time
    % PLV = abs(mean(exp(1i * (phi_it - phi_pfc)), across trials))
    PLV_it_pfc(ci,:,:) = squeeze(abs(nanmean(exp(1i*(phi_it - phi_pfc)), 1)));  % dims: freq x time

    % Save phase difference (radians) per trial (trial x freq x time)
    phaseDiff{ci} = (phi_it - phi_pfc);

    % Shuffling: random trial pairing to estimate baseline PLV
    nTrialsCond = size(phi_it,1);
    parfor itri = 1:N_shuff   % use parfor if Parallel Toolbox available
        % Draw random permutations of trial indices for each side
        ind_sh1 = randperm(nTrialsCond);
        ind_sh2 = randperm(nTrialsCond);
        % compute shuffled PLV (mean across trials after random pairing)
        tmp = abs(nanmean(exp(1i*(phi_it(ind_sh1,:,:) - phi_pfc(ind_sh2,:,:))), 1));
        PLV_pfc_it_shuff(ci, itri, :, :) = squeeze(tmp);   % freq x time
    end

    % Compute baseline-corrected PLV: observed minus mean(shuffled)
    mean_shuff = squeeze(nanmean(PLV_pfc_it_shuff(ci,:,:,:), 2));  % freq x time
    PLV_shuf_c(ci,:,:) = squeeze(PLV_it_pfc(ci,:,:)) - mean_shuff;
end

% Save results (same variables you used)
% savePath = fullfile('G:\Neural_DATA\LFP\PLV\Data', ['_Shufl_onTime_' R.name '.mat']);
% save(savePath, "PLV_shuf_c", "phaseDiff", "PLV_it_pfc", "PLV_pfc_it_shuff", "sig_pfc_wav", "sig_it_wav", "-v7.3");
% ---- Plotting: Min–Max Normalized TF PLV images + band timecourses ----
clc;

% Define canonical frequency bands
bands = struct('name', {'Theta','Alpha','Beta','LowGamma','HighGamma'}, ...
               'range', {[4 8], [8 13], [13 30], [31 45], [55 120]});

timeVec = linspace(-200,700,900);    % time vector (ms)

% ---- Apply Min–Max normalization per frequency ----
PLV_norm = zeros(size(PLV_shuf_c));   % same shape: cond x freq x time

for ci = 1:nConds
    for f = 1:length(freq)
        rawVals = squeeze(PLV_shuf_c(ci,f,:));
        minVal = min(rawVals, [], 'omitnan');
        maxVal = max(rawVals, [], 'omitnan');

        if maxVal > minVal  % avoid divide-by-zero
            PLV_norm(ci,f,:) = (rawVals - minVal) / (maxVal - minVal);
        else
            PLV_norm(ci,f,:) = zeros(size(rawVals)); % flat if no variation
        end
    end
end

% ---- Plot Smoothed Normalized TF PLV Maps ----
figure('Position',[100 100 1200 800]);

% Define smoothing window (adjust for desired smoothness)
smoothWin_time = 20;  % in samples along time
smoothWin_freq = 10;  % in samples along frequency
h = fspecial('gaussian', [smoothWin_freq smoothWin_time], 5); % Gaussian filter

for ci = 1:nConds
    % Extract matrix for this condition
    PLVmat = squeeze(PLV_norm(ci,:,:));
    
    % Apply 2D smoothing
    PLV_smooth = imfilter(PLVmat, h, 'replicate');
    
    % Plot smoothed heatmap
    subplot(2,2,ci);
    imagesc(timeVec, freq, PLV_smooth);
    axis xy;
    xlabel('Time (ms)');
    ylabel('Frequency (Hz)');
    ylim([0 20]);
    xlim([-100 600]);
    title(sprintf('Condition %d: Smoothed Normalized PLV', ci));
    colorbar;
    caxis([0 1]); % normalized scale
end


% ---- Plot band-specific normalized PLV over time ----
figure('Position',[100 100 1200 600]);
bandPLV_all = struct();  % store for saving later

for ci = 1:nConds
    subplot(nConds,1,ci);
    hold on;
    for b = 1:length(bands)
        fmask = freq >= bands(b).range(1) & freq <= bands(b).range(2);
        if ~any(fmask)
            continue;
        end
        bandPLV = squeeze(PLV_norm(ci, fmask, :));   % freq_in_band x time
        bandMean = nanmean(bandPLV, 1);              % mean across frequencies
        plot(timeVec, movmean(bandMean,20), 'LineWidth', 1.5);

        % --- Store mean PLV per band and condition for later use ---
        bandPLV_all(ci).(bands(b).name) = bandMean;
    end
    xlabel('Time (ms)');
    ylabel('Normalized PLV (0–1)');
    title(sprintf('Condition %d: Normalized PLV by Band', ci));
    legend({bands.name}, 'Location','northeastoutside');
    grid on;
    ylim([0 1]);
    xlim([-100 600]);
end

% ---- Extract and save Beta band PLV ----
% Get Beta band PLV (13–30 Hz) for all conditions
Beta_PLV_all = cell(nConds,1);
for ci = 1:nConds
    if isfield(bandPLV_all(ci), 'Beta')
        Beta_PLV_all{ci} = bandPLV_all(ci).Beta;
    else
        Beta_PLV_all{ci} = [];
    end
end

% Save to .mat file
save('Beta_PLV_across_time_late_high_FR.mat', 'Beta_PLV_all', 'timeVec', 'bands');

disp('✅ Beta PLV timecourses saved to "Beta_PLV_across_time_late_high_FR.mat"');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---- Optional: Circular mean of phase difference in a time-frequency window ----
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example: compute mean angle in a time window (e.g., 150:200 ms relative) and beta band
time_window = find(timeVec>=100 & timeVec<=200);   % adjust to desired ms window
band_idx = freq >= 13 & freq <= 30;                % beta
figure;
for ci = 1:nConds
    % phaseDiff{ci} is trial x freq x time
    pd = phaseDiff{ci}(:, band_idx, time_window);   % trials x freq_in_band x time_in_window
    % collapse freq & time, then compute circular mean
    pd_vec = pd(:);
    mean_ang = angle(nanmean(exp(1i*pd_vec)));
    % Plot as polar
    subplot(2,2,ci);
    polarplot([mean_ang mean_ang], [0 1], 'LineWidth', 3);
    title(sprintf('Cond %d: Mean phase diff (beta, %d-%d ms)', ci, timeVec(time_window(1)), timeVec(time_window(end))));
end

% Done
fprintf('PLV computation + plots complete. Results saved to:\n%s\n', savePath);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---- Load Data ---- early 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataLow_early  = load('Beta_PLV_across_time_early_low.mat');
dataHigh_early = load('Beta_PLV_across_time_early_high.mat');
dataLow_late  = load('Beta_PLV_across_time_late_low.mat');
dataHigh_late = load('Beta_PLV_across_time_late_high.mat');

Beta_PLV_low_all  = dataLow.Beta_PLV_all;
Beta_PLV_high_all = dataHigh.Beta_PLV_all;

timeVec = linspace(-200,700,900);

% ---- Parameters ----
baselineWindow = [-200 0];      % baseline period
timeWindowEarly = [0 300];
timeWindowLate  = [300 600];

% Indices
tmaskBase  = timeVec >= baselineWindow(1) & timeVec <= baselineWindow(2);
tmaskEarly = timeVec >= timeWindowEarly(1) & timeVec <= timeWindowEarly(2);
tmaskLate  = timeVec >= timeWindowLate(1)  & timeVec <= timeWindowLate(2);

% ---- Baseline Correction ----
for n = 1:numel(Beta_PLV_low_all)
    baseMean_low  = mean(Beta_PLV_low_all{n}(tmaskBase), 'omitnan');
    baseMean_high = mean(Beta_PLV_high_all{n}(tmaskBase), 'omitnan');
    Beta_PLV_low_all{n}  = Beta_PLV_low_all{n}  - baseMean_low;
    Beta_PLV_high_all{n} = Beta_PLV_high_all{n} - baseMean_high;
end

% ---- Compute mean across neurons for timecourse ----
Beta_PLV_mean_low  = mean(cell2mat(cellfun(@(x) x(:), Beta_PLV_low_all, 'UniformOutput', false)), 2);
Beta_PLV_mean_high = mean(cell2mat(cellfun(@(x) x(:), Beta_PLV_high_all, 'UniformOutput', false)), 2);

% ---- Plot Beta PLV Across Time (Line Plot) ----
figure('Position',[100 100 1200 500]);
plot(timeVec, movmean(Beta_PLV_low_all{n},100), 'b', 'LineWidth', 1.5); hold on;
plot(timeVec, movmean(Beta_PLV_high_all{n},100), 'r', 'LineWidth', 1.5);
xline(0, '--k'); xline(300, '--k'); xline(600, '--k');
grid on;
xlabel('Time (ms)');
ylabel('Baseline-Corrected Beta PLV');
legend({'Low Neural Var','High Neural Var'}, 'Location','best');
title('Beta-band PLV Across Time (Smoothed)');

% ---- Compute mean & SE for Bar Plot ----
Conditions = numel(Beta_PLV_low_all);
meanLow_early  = zeros(Conditions,1);
meanLow_late   = zeros(Conditions,1);
meanHigh_early = zeros(Conditions,1);
meanHigh_late  = zeros(Conditions,1);

for n = 1:Conditions
    meanLow_early(n)  = mean(Beta_PLV_low_all{n}(tmaskEarly), 'omitnan');
    meanLow_late(n)   = mean(Beta_PLV_low_all{n}(tmaskLate), 'omitnan');
    meanHigh_early(n) = mean(Beta_PLV_high_all{n}(tmaskEarly), 'omitnan');
    meanHigh_late(n)  = mean(Beta_PLV_high_all{n}(tmaskLate), 'omitnan');
end

means = [
    mean(meanLow_early),  mean(meanHigh_early);
    mean(meanLow_late),   mean(meanHigh_late)
];
SEs = [
    std(meanLow_early)/sqrt(Conditions),  std(meanHigh_early)/sqrt(Conditions);
    std(meanLow_late)/sqrt(Conditions),   std(meanHigh_late)/sqrt(Conditions)
];

% ---- Wilcoxon rank-sum tests ----
[pEarly,~] = ranksum(meanLow_early, meanHigh_early);
[pLate,~]  = ranksum(meanLow_late,  meanHigh_late);

% ---- Plot Bar Plot with Error Bars and Significance ----
figure('Position',[100 100 800 500]); hold on;
barHandles = bar(means, 'grouped');
barHandles(1).FaceColor = [0.3 0.5 1];   % LowVar
barHandles(2).FaceColor = [1 0.4 0.4];   % HighVar

ngroups = size(means,1);
nbars = size(means,2);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, means(:,i), SEs(:,i), 'k', 'linestyle','none', 'LineWidth',1.2);
end

set(gca, 'XTickLabel', {'Early (0–300 ms)', 'Late (300–600 ms)'});
xlabel('Time Window'); ylabel('Mean Beta PLV');
legend({'Low Neural Var','High Neural Var'},'Location','northwest');
title('Beta-band PLV: Low vs High Neural Variability');
grid on;

% Significance stars
for w = 1:2
    if w==1, pval = pEarly; y = max(means(1,:) + SEs(1,:)) + 0.05; x = 1; end
    if w==2, pval = pLate;  y = max(means(2,:) + SEs(2,:)) + 0.05; x = 2; end
    if pval < 0.001, stars='***';
    elseif pval < 0.01, stars='**';
    elseif pval < 0.05, stars='*';
    else stars='n.s.'; end
    text(x, y, stars, 'HorizontalAlignment','center', 'FontSize',14);
end

ylim([min(means(:)-SEs(:))-0.05, max(means(:)+SEs(:))+0.15]);

% ---- Save Results ----
save('Mean_Beta_PLV_withBaseline_Wilcoxon.mat', 'means', 'SEs', 'pEarly', 'pLate', ...
    'baselineWindow', 'timeWindowEarly', 'timeWindowLate', 'Beta_PLV_mean_low','Beta_PLV_mean_high');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---- Load Data for comparision Beta band PLV ----
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataLow_early  = load('Beta_PLV_across_time_early_low_FR.mat');
dataHigh_early = load('Beta_PLV_across_time_early_high_FR.mat');
dataLow_late   = load('Beta_PLV_across_time_late_low_FR.mat');
dataHigh_late  = load('Beta_PLV_across_time_late_high_FR.mat');

Beta_low_early  = dataLow_early.Beta_PLV_all;
Beta_high_early = dataHigh_early.Beta_PLV_all;
Beta_low_late   = dataLow_late.Beta_PLV_all;
Beta_high_late  = dataHigh_late.Beta_PLV_all;

timeVec = linspace(-200,700,900);

% ---- Define Time Windows ----
timeWindowEarly = [60 290];
timeWindowLate  = [300 500];

tmaskEarly = timeVec >= timeWindowEarly(1) & timeVec <= timeWindowEarly(2);
tmaskLate  = timeVec >= timeWindowLate(1)  & timeVec <= timeWindowLate(2);

% ---- Compute mean ± SE across neurons (no baseline correction) ----
Conditions = numel(Beta_low_early);

meanLow_early  = zeros(Conditions,1);
meanHigh_early = zeros(Conditions,1);
meanLow_late   = zeros(Conditions,1);
meanHigh_late  = zeros(Conditions,1);

for n = 1:Conditions
    meanLow_early(n)  = mean(Beta_low_early{n}(tmaskEarly), 'omitnan');
    meanHigh_early(n) = mean(Beta_high_early{n}(tmaskEarly), 'omitnan');
    meanLow_late(n)   = mean(Beta_low_late{n}(tmaskLate), 'omitnan');
    meanHigh_late(n)  = mean(Beta_high_late{n}(tmaskLate), 'omitnan');
end

means = [
    mean(meanLow_early),  mean(meanHigh_early);
    mean(meanLow_late),   mean(meanHigh_late)
];

SEs = [
    std(meanLow_early)/sqrt(Conditions),  std(meanHigh_early)/sqrt(Conditions);
    std(meanLow_late)/sqrt(Conditions),   std(meanHigh_late)/sqrt(Conditions)
];

% ---- Wilcoxon rank-sum tests ----
[pEarly,~] = ranksum(meanLow_early, meanHigh_early);
[pLate,~]  = ranksum(meanLow_late,  meanHigh_late);

% ---- Plot Bar Graph with Error Bars ----
figure('Position',[100 100 800 500]); hold on;
barHandles = bar(means,'grouped');
barHandles(1).FaceColor = [0.3 0.5 1];   % LowVar
barHandles(2).FaceColor = [1 0.4 0.4];   % HighVar

ngroups = size(means,1);
nbars = size(means,2);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth/(2*nbars);
    errorbar(x, means(:,i), SEs(:,i), 'k','linestyle','none','LineWidth',1.2);
end

% Hide x-tick labels but keep time frame text above bars
set(gca,'XTickLabel',[]);
ylim([min(means(:)-SEs(:))-0.05, max(means(:)+SEs(:))+0.15]);

% Add time frame text above the groups
text(1, max(means(1,:) + SEs(1,:)) + 0.08, 'Early (60–290 ms)', 'HorizontalAlignment','center', 'FontSize',12);
text(2, max(means(2,:) + SEs(2,:)) + 0.08, 'Late (300–500 ms)',  'HorizontalAlignment','center', 'FontSize',12);

% Add significance stars
for w = 1:2
    if w==1, pval = pEarly; y = max(means(1,:) + SEs(1,:)) + 0.05; x = 1; end
    if w==2, pval = pLate;  y = max(means(2,:) + SEs(2,:)) + 0.05; x = 2; end
    if pval < 0.001, stars='***';
    elseif pval < 0.01, stars='**';
    elseif pval < 0.05, stars='*';
    else stars='n.s.'; end
    text(x, y, stars,'HorizontalAlignment','center','FontSize',14);
end

xlabel('Time Window'); ylabel('Mean Beta PLV');
legend({'Low Firing rate','High Firing rate'},'Location','northwest');
title('Beta-band PLV: Low vs High Firing Rate');

% ---- Save results ----
save('Mean_Beta_PLV_EarlyLate_Wilcoxon_noBaseline.mat','means','SEs','pEarly','pLate',...
    'timeWindowEarly','timeWindowLate');



% ============================================
%% ===== Robust loop to compute Beta PLV per session =====
freq = [1:1:45 55:1:160];
fs = 1000; pad = 1;
N_shuff = 500;
bands = struct('name', {'Theta','Alpha','Beta','LowGamma','HighGamma'}, ...
               'range', {[4 8], [8 13], [13 30], [31 45], [55 120]});
timeVec = linspace(-200,700,900);

for ss = 1:110
    fprintf('Processing session (neuron) %d / %d\n', ss, length(SpikeTrain_it_all));

    sig_it  = SpikeTrain_it_all(ss).data_it;     
    sig_pfc = SpikeTrain_it_all(ss).data_pfc;
    if isempty(sig_it) || isempty(sig_pfc)
        warning('Session %d: missing data_it or data_pfc — skipping.', ss);
        continue;
    end

    artif_it  = lfp_artifact_remove(sig_it, 0);
    artif_pfc = lfp_artifact_remove(sig_pfc, 0);
    artif = artif_it | artif_pfc;

    [sig_it_an, fr]   = ndass_wavelet_np(sig_it, freq, fs, pad);
    [sig_pfc_an, fr]  = ndass_wavelet_np(sig_pfc, freq, fs, pad);

    if isempty(sig_it_an) || isempty(sig_pfc_an)
        warning('Session %d: invalid wavelet output — skipping.', ss);
        continue;
    end

    nTrials_total = size(sig_it_an,1);
    nFreq = size(sig_it_an,2);
    nTime = size(sig_it_an,3);

    ind_h = highVarTrials_late{ss};
    if isempty(ind_h)
        warning('Session %d: no valid trials — skipping.', ss);
        continue;
    end
    if islogical(ind_h) && length(ind_h) ~= nTrials_total
        ind_h = find(ind_h);
    end

    phi_it  = angle(sig_it_an(ind_h,:,:));
    phi_pfc = angle(sig_pfc_an(ind_h,:,:));

    if isempty(phi_it) || isempty(phi_pfc)
        warning('Session %d: no valid phases — skipping.', ss);
        continue;
    end

    obs = nanmean(exp(1i*(phi_it - phi_pfc)),1);
    PLV_it_pfc = squeeze(abs(obs));

    nTrialsCond = size(phi_it,1);
    if nTrialsCond < 2
        warning('Session %d: not enough trials — skipping.', ss);
        continue;
    end

    PLV_shuff = nan(N_shuff, nFreq, nTime);
    parfor itri = 1:N_shuff
        ind_sh1 = randperm(nTrialsCond);
        ind_sh2 = randperm(nTrialsCond);
        tmp = nanmean(exp(1i*(phi_it(ind_sh1,:,:) - phi_pfc(ind_sh2,:,:))), 1);
        PLV_shuff(itri,:,:) = squeeze(abs(tmp));
    end
    mean_shuff = squeeze(nanmean(PLV_shuff,1));
    PLV_shuf_c = PLV_it_pfc - mean_shuff;

    PLV_norm = zeros(size(PLV_shuf_c));
    for f = 1:nFreq
        rawVals = squeeze(PLV_shuf_c(f,:));
        minVal = min(rawVals,[],'omitnan');
        maxVal = max(rawVals,[],'omitnan');
        if maxVal > minVal
            PLV_norm(f,:) = (rawVals - minVal) / (maxVal - minVal);
        else
            PLV_norm(f,:) = zeros(1,nTime);
        end
    end

    fmask_beta = freq >= 13 & freq <= 30;
    BetaPLV_vec = nanmean(PLV_norm(fmask_beta,:),1);

    % ✅ safe: modify structure in serial loop
    SpikeTrain_it_all(ss).BetaPLV_HighVar_Late = BetaPLV_vec;
end

disp('✅ Beta band PLV vectors computed and saved in SpikeTrain_it_all');

%%

save('PLV_Beta_HighLOwVar_early.mat','SpikeTrain_it_all');
%% ============================================
%   SCATTER PLOTS: Session-level PLV (Early & Late)
% ============================================

figure('Position',[200 100 1000 450]);

% ---- Early Time Window Scatter ----
subplot(1,2,1);
scatter(meanHigh_early, meanLow_early, 80, 'filled', ...
        'MarkerFaceColor',[0.2 0.6 0.9], 'MarkerEdgeColor','k');
hold on;
plot([min([meanHigh_early; meanLow_early]) max([meanHigh_early; meanLow_early])], ...
     [min([meanHigh_early; meanLow_early]) max([meanHigh_early; meanLow_early])], ...
     'r--', 'LineWidth',1.5);
xlabel('High Variability (Mean PLV)');
ylabel('Low Variability (Mean PLV)');
title('Early Time Window (60–290 ms)');
axis equal;
grid on;
legend('Sessions','y = x','Location','best');
xlim([min([meanLow_early; meanHigh_early]) max([meanLow_early; meanHigh_early])]);
ylim([min([meanLow_early; meanHigh_early]) max([meanLow_early; meanHigh_early])]);

% ---- Late Time Window Scatter ----
subplot(1,2,2);
scatter(meanHigh_late, meanLow_late, 80, 'filled', ...
        'MarkerFaceColor',[0.9 0.5 0.2], 'MarkerEdgeColor','k');
hold on;
plot([min([meanHigh_late; meanLow_late]) max([meanHigh_late; meanLow_late])], ...
     [min([meanHigh_late; meanLow_late]) max([meanHigh_late; meanLow_late])], ...
     'r--', 'LineWidth',1.5);
xlabel('High Variability (Mean PLV)');
ylabel('Low Variability (Mean PLV)');
title('Late Time Window (300–500 ms)');
axis equal;
grid on;
legend('Sessions','y = x','Location','best');
xlim([min([meanLow_late; meanHigh_late]) max([meanLow_late; meanHigh_late])]);
ylim([min([meanLow_late; meanHigh_late]) max([meanLow_late; meanHigh_late])]);

sgtitle('Session-level Beta-band PLV: High vs Low Variability');
nEarly_higherHigh = sum(meanHigh_early > meanLow_early);
nLate_higherHigh  = sum(meanHigh_late > meanLow_late);

fprintf('Early window: %d/%d sessions had higher PLV in high variability\n', ...
    nEarly_higherHigh, numel(meanHigh_early));
fprintf('Late window: %d/%d sessions had higher PLV in high variability\n', ...
    nLate_higherHigh, numel(meanHigh_late));


%% ==========================
%  SCATTER PLOT (Low vs High Variability)
% ==========================

% ---- Early ----
figure('Position',[200 200 500 450]); hold on;

% Scatter points
scatter(meanHigh_early, meanLow_early, 60, 'k', 'filled', 'MarkerFaceAlpha',0.7);
plot([min([meanHigh_early; meanLow_early]) max([meanHigh_early; meanLow_early])], ...
     [min([meanHigh_early; meanLow_early]) max([meanHigh_early; meanLow_early])], ...
     'k-', 'LineWidth',1); % equality line

% Axis labels
xlabel('High Variability PLV');
ylabel('Low Variability PLV');
title('Early Time Window (60–290 ms)');
axis square; box on;

% Count how many points above/below y=x
above = sum(meanLow_early > meanHigh_early);
below = sum(meanHigh_early > meanLow_early);
text(min(xlim)+0.05*range(xlim), max(ylim)-0.1*range(ylim), sprintf('N=%d', above), 'FontSize',12);
text(max(xlim)-0.3*range(xlim), min(ylim)+0.05*range(ylim), sprintf('N=%d', below), 'FontSize',12, 'Color','r');

% Wilcoxon test
[pEarly,~] = ranksum(meanLow_early, meanHigh_early);
text(max(xlim)-0.5*range(xlim), min(ylim)+0.15*range(ylim), sprintf('p = %.3f', pEarly), 'FontSize',12);

% ---- Small inset histogram of PLV difference ----
axes('Position',[0.65 0.65 0.25 0.25]); % inset position
diffEarly = meanHigh_early - meanLow_early;
histogram(diffEarly, 'FaceColor',[0.3 0.3 0.3], 'EdgeColor','none');
hold on;
plot([0 0], ylim, 'k-', 'LineWidth',1);
xlabel('\DeltaPLV');
set(gca, 'FontSize',8, 'Box','off');

sgtitle('Beta-band PLV: Early window');


% ---- Late ----
figure('Position',[750 200 500 450]); hold on;

scatter(meanHigh_late, meanLow_late, 60, 'k', 'filled', 'MarkerFaceAlpha',0.7);
plot([min([meanHigh_late; meanLow_late]) max([meanHigh_late; meanLow_late])], ...
     [min([meanHigh_late; meanLow_late]) max([meanHigh_late; meanLow_late])], ...
     'k-', 'LineWidth',1);

xlabel('High Variability PLV');
ylabel('Low Variability PLV');
title('Late Time Window (300–500 ms)');
axis square; box on;

above = sum(meanLow_late > meanHigh_late);
below = sum(meanHigh_late > meanLow_late);
text(min(xlim)+0.05*range(xlim), max(ylim)-0.1*range(ylim), sprintf('N=%d', above), 'FontSize',12);
text(max(xlim)-0.3*range(xlim), min(ylim)+0.05*range(ylim), sprintf('N=%d', below), 'FontSize',12, 'Color','r');

[pLate,~] = ranksum(meanLow_late, meanHigh_late);
text(max(xlim)-0.5*range(xlim), min(ylim)+0.15*range(ylim), sprintf('p = %.3f', pLate), 'FontSize',12);

axes('Position',[0.65 0.65 0.25 0.25]); % inset
diffLate = meanHigh_late - meanLow_late;
histogram(diffLate, 'FaceColor',[0.3 0.3 0.3], 'EdgeColor','none');
hold on;
plot([0 0], ylim, 'k-', 'LineWidth',1);
xlabel('\DeltaPLV');
set(gca, 'FontSize',8, 'Box','off');

sgtitle('Beta-band PLV: Late window');
