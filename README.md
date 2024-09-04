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
