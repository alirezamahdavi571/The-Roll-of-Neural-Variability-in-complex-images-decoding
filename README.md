![FIGURE1](https://github.com/user-attachments/assets/dea2a0c6-e449-43cb-bf4c-47bcd2a5a501)
# Stimuli and task
The task paradigm required maintenance of eye fixation within a 2° radius around the center of a monitor (Asus VG248QE: 24in, resolution 1920x1080, 144 Hz refresh rate). Upon ensuring the eye fixation, a stimulus is displayed in place of the fixation point for a 50-ms duration over a gray background, followed by a 600-ms blank interval. The stimuli were scaled to fit in a 7° by 7° square at the center of the monitor. Monkeys got rewarded for holding fixation every 1.5 seconds with droppings of sweet juice. The stimulus set contains 155 pictures from 4 different categories: face, body, artificial, and natural. 
Behavioral and recording method
We conducted experiments using acute extracellular recording in behaving rhesus macaques. Using two single tungsten electrodes (A-M Systems, FHC) simultaneously, we recorded the broadband neural activity of the prefrontal and inferior temporal cortices. Since the receptive field (RF) of IT neurons is relatively large and the RF of PFC neurons is not well-defined, we determined the spatial preference based on the PFC neuron’s RF. Before starting the main task, the RF of PFC neurons was determined by performing an RF-mapping task and monitoring the responses of the neurons to the spatial location of the stimulus. A reference grid, with holes 1 mm apart, was used inside all the recording chambers to guide electrode penetrations and localize them relative to structural MRI images. The neural response was amplified, quantized, and saved using a 40KHz data acquisition system (NikTek Systems).
# Neural data preprocessing
 Wide-band neural data was initially band-pass filtered between 300 to 3000 Hz to extract the high-frequency signals. By setting a multiplier of the median of the resulting waveform, the spikes were detected using the ROSS Toolbox 28. Using different electrodes for recording from different monkeys, we did not use the sorted spikes in the analysis and used the multiunit activity instead. Followed by spike detection, we chose a subset of recorded neurons using a selectivity criteria by which a unit was discarded if its post-stimulus average firing rate didn’t change significantly compared to its baseline activity for all of the four categories. The statistical significance was measured by a Wilcoxon signed rank test with a 5% confidence level.

#### Files Structure

```bash
SRC
    ├── Fano factor
                    ├── Category based
                    └── Stimuli based
    ├── GLM
    ├── Information theory
                        ├── Mutual Information
                        └── Factor Analysis
    ├── PSTH
    └── Machine learning methods
                            ├── LDA
                            ├── Generalization (Time-Time decoding)
                            └── SVM                          -       
```
The diagram below shows the entire project and the data flow:
![jhidnjsbdw copy](https://github.com/user-attachments/assets/16d9f453-3d1c-40dd-b19f-e04745ab9077)
## How to use codes in this Project
There are some exmple here that show how to use codes and functions in matlab.
### Load data 
I have 859 neurons that any neuron contains a raster matrix from spiking data`SpikeTrain_it_all`. Please consider this note that the input data is sorted.
for first step we load data in matlab:

```matlab
% List all folders in the current directory
folderList = dir('*.*'); % List all files and folders
folderList = folderList([folderList.isdir]); % Filter out only folders
folderNames = {folderList.name}; % Extract folder names
disp('Folders in the current directory:');
disp(folderNames);

% Initialize structure array to store data from each folder
SpikeTrain_it_all = struct('data', []);

% Open each folder
for i = 1:length(folderNames)
    if ~strcmp(folderNames{i}, '.') && ~strcmp(folderNames{i}, '..')
        folderPath = fullfile(pwd, folderNames{i} , '\Trial'); % Get full path of the folder
        disp(['Opening folder: ', folderPath]);
        cd(folderPath); % Change current directory to the folder
        % Load the mu_it.mat file from the Trial folder
        mu_it_data = load('mu_it.mat');
        cm = load('cm.mat');
        % Store the loaded data into the structure array
        SpikeTrain_it_all(i).data = mu_it_data;
        SpikeTrain_it_all(i).cm = cm;
        % Return to the parent directory
        cd('../..'); 
    end
end                       -       
```
### seperate face data and non face cm and trials conditions
for have `SVM` classifier on data raster for any classes should be seperated.
```matlab
%% seperate face data and non face cm and trials conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
SpikeTrain_it_all = seperateStimulus(SpikeTrain_it_all,face_labels,'face',number_of_neurons);
SpikeTrain_it_all = seperateStimulus(SpikeTrain_it_all,body_labels,'body',number_of_neurons);
SpikeTrain_it_all = seperateStimulus(SpikeTrain_it_all,natural_labels,'natural',number_of_neurons);
SpikeTrain_it_all = seperateStimulus(SpikeTrain_it_all,artifact_labels,'artifact',number_of_neurons);
SpikeTrain_it_all = seperateStimulus(SpikeTrain_it_all,nonface_labels,'nonface',number_of_neurons);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
```
#### calculate PSTH(firing rate)
features should be extracted from neuron's firing rate (PSTH) 
```matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% default is seperation between face and non-face conditions
% first step is seperate raster for any stimulus
% second step is PSTH that has shpe 
%                                    = (number_of_neurons ,number_of_labels_class,time linspace)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
SpikeTrain_it_all = PSTH(SpikeTrain_it_all,number_of_neurons,'face',sliding_step,window_length,number_of_time_slices);
SpikeTrain_it_all = PSTH(SpikeTrain_it_all,number_of_neurons,'nonface',sliding_step,window_length,number_of_time_slices);
SpikeTrain_it_all = PSTH(SpikeTrain_it_all,number_of_neurons,'body',sliding_step,window_length);
SpikeTrain_it_all = PSTH(SpikeTrain_it_all,number_of_neurons,'artifact',sliding_step,window_length);
SpikeTrain_it_all = PSTH(SpikeTrain_it_all,number_of_neurons,'natural',sliding_step,window_length);

```
### SVM
``` matlab
Face_Mean_PSTH_data = zeros(number_of_neurons,min_stimulus,number_of_time_slices);
body_Mean_PSTH_data = zeros(number_of_neurons,min_stimulus,number_of_time_slices);

for i = 3:number_of_neurons
    Face_Mean_PSTH_data(i,:,:) = SpikeTrain_it_all(i).Face_Mean_PSTH_data(1:min_stimulus,:);
    body_Mean_PSTH_data(i,:,:) = SpikeTrain_it_all(i).NonFace_Mean_PSTH_data(1:min_stimulus,:);
end

save('Face_Mean_PSTH_data.mat' , 'Face_Mean_PSTH_data' ,'-v7.3' );
save('nonFace_Mean_PSTH_data.mat' , 'body_Mean_PSTH_data' ,'-v7.3' );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SVM classifier
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
for ts = 1:number_of_time_slices
    face_features = transpose(Face_Mean_PSTH_data(:, :, ts));
    nonface_features = transpose(body_Mean_PSTH_data(:, :, ts));

    % Create labels
    face_labels = -1 * ones(size(face_features, 1), 1);  % -1 for face
    nonface_labels = ones(size(nonface_features, 1), 1);     % +1 for non-faceanydesk

    % Stack features and labels
    features = [face_features; nonface_features];
    labels = [face_labels; nonface_labels];
    disp(['Accuracy for time slice ' num2str(ts)]);
    out(ts).out = gen_fx_get_svm(labels,features,0.8,iteration);
end


for ts = 1:number_of_time_slices
    accuracy_w(ts,1) = mean(out(ts).out.pt,1)* 100;
end
```
accuracy across time slices can be watched:
![image](https://github.com/user-attachments/assets/2a8e8a49-d625-4a69-b7f6-58fadb034815)

### Fano factor category based
# Fano Factor Calculation

The **Fano Factor** is a measure of the dispersion of a probability distribution and is defined as the ratio of the variance to the mean of a dataset. It is commonly used in fields such as neuroscience, physics, and statistics to assess the variability of spike counts or other data.

#### Formula

The Fano Factor (FF) is given by:

$$
\text{Fano Factor (FF)} = \frac{\sigma^2}{\mu}
$$

where:
- $$(\sigma^2\)$$ is the variance of the data
- $$(\mu\)$$ is the mean of the data

#### Usage

To calculate the Fano Factor:
1. Collect the dataset for which you want to calculate the Fano Factor.
2. Compute the mean $$(\(\mu\))$$ of the dataset.
3. Compute the variance $$(\(\sigma^2\))$$ of the dataset.
4. Use the formula above to compute the Fano Factor.

#### Example

Suppose we have a dataset representing spike counts in a neuroscience experiment:
```matlab
[mean_vec_face,var_vec_face,SpikeTrain_it_all] = catagoryBasedFano(SpikeTrain_it_all,'face',face_labels,mean_vec_face,var_vec_face,...
    number_of_neurons,sliding_step,window_length,number_of_time_slices);

[mean_vec_body,var_vec_body,SpikeTrain_it_all] = catagoryBasedFano(SpikeTrain_it_all,'body',body_labels,mean_vec_body,var_vec_body,...
    number_of_neurons,sliding_step,window_length,number_of_time_slices);

[mean_vec_natural,var_vec_natural,SpikeTrain_it_all] = catagoryBasedFano(SpikeTrain_it_all,'natural',natural_labels,mean_vec_natural,var_vec_natural,...
    number_of_neurons,sliding_step,window_length,number_of_time_slices);

[mean_vec_artifact,var_vec_artifact,SpikeTrain_it_all] = catagoryBasedFano(SpikeTrain_it_all,'artifact',artifact_labels,mean_vec_artifact,var_vec_artifact,...
    number_of_neurons,sliding_step,window_length,number_of_time_slices);

[mean_vec_nonface,var_vec_nonface,SpikeTrain_it_all] = catagoryBasedFano(SpikeTrain_it_all,'nonface',nonface_labels,mean_vec_nonface,var_vec_nonface,...
    number_of_neurons,sliding_step,window_length,number_of_time_slices);

for i = 1:number_of_time_slices
    [fanofactor_face(i),~,~] = regression(mean_vec_face(:,i), var_vec_face(:,i),'one');
    [fanofactor_artfact(i),~,~] = regression(mean_vec_artifact(:,i), var_vec_artifact(:,i),'one');
    
    [fanofactor_body(i),~,~] = regression(mean_vec_body(:,i), var_vec_body(:,i),'one');
    [fanofactor_natural(i),~,~] = regression(mean_vec_natural(:,i), var_vec_natural(:,i),'one');
    [fanofactor_nonface(i),~,~] = regression(mean_vec_nonface(:,i), var_vec_nonface(:,i),'one');
end
```

please refer to Churchland's paper to inform about mean match fanofactor: `https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2828350/`

### Fano factor stimuli based

```matlab
All_stimuli = [artifact_labels ; body_labels ; natural_labels ; face_labels];
indexes = [];
mean_vec = zeros(number_of_neurons , length(All_stimuli) , number_of_time_slices);
var_vec = zeros(number_of_neurons , length(All_stimuli) , number_of_time_slices);
[mean_vec,var_vec] = ffstmbased(SpikeTrain_it_all , mean_vec,var_vec,All_stimuli,number_of_neurons,number_of_time_slices,sliding_step,window_length);
```

### Entropy
we expect that Entropy behiaves like FF:
``` matlab
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate information 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

entropy_nonface = entroInf(mean_vec_nonface,var_vec_nonface);
entropy_face = entroInf(mean_vec_face,var_vec_face);
time = linspace(-200,700,number_of_time_slices);

figure()
subplot(2,1,1)
p1 = plot(time , entropy_face, 'LineWidth',5);
grid on
hold on
title('entropy of Fano factor');ylabel('ratio');xlabel('time');
p2 = plot(time , entropy_nonface, 'LineWidth',5);
legend([p1,p2] , {'face' , 'nonface'} , 'FontSize',12)

subplot(2,1,2)
p1 = plot(time , fanofactor_face, 'LineWidth',5);
grid on
hold on
title('Fanofactor');ylabel('ratio');xlabel('time');
p2 = plot(time , fanofactor_nonface, 'LineWidth',5);
legend([p1,p2] , {'face' , 'nonface'} , 'FontSize',12)
```
Entropy for 2 classes can be looked:
![image](https://github.com/user-attachments/assets/8d6ca83d-a9db-45cb-a4e9-37feed8bd865)










