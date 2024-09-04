# Spiking-codes
![2](https://github.com/user-attachments/assets/c755ad02-fb95-4f5e-882b-2064fa8efb3c)


In this project, we have put almost all the codes needed for spiking data processing, which include the following sections:

#### Files Structure

```bash
SRC
    ├── Fano factor
    ├── GLM
    ├── Information theory
    ├── PSTH
    └── Machine learning methods
                            └──SVM                          -       
```
The diagram below shows the entire project and the data flow:
![jhidnjsbdw copy](https://github.com/user-attachments/assets/16d9f453-3d1c-40dd-b19f-e04745ab9077)
## How to use codes in this Project
There are some exmple here that show how to use codes and functions in matlab.
#### Load data 
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
#### seperate face data and non face cm and trials conditions
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
#### SVM
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











