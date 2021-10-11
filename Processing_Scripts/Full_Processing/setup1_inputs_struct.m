% author: Nicole Bienert, Rohan Sanda
% date: 7/12/21
% Purpose: To record the initial variable settings into a struct object so
% that they can be automatically pulled in each script. This way, we can
% ensure the continuity of numerical constants throughout all phases of the
% processing chain.

clc
close all
clear all

%Get user inputs
prompt1 = sprintf('Would you like to process data from Greenland or Antarctica? \nEnter 1 or 2 to select your choice \n   1: Greenland \n   2: Antarctica\n'); 
user_G_or_A = input(prompt1);
if(user_G_or_A ~= 1 && user_G_or_A ~= 2)
    error('User input not valid.')
end

if (user_G_or_A ==2)
    prompt2 = sprintf('\nWhich Antenna Separation? \nEnter 1, 2, 3, or 4 to select your choice. \n   1: 300 m \n   2: 600 m \n   3: 900 m \n   4: 1300 m\n');
    user_offset = input(prompt2);
    if(user_offset ~= 1 && user_offset ~= 2 && user_offset~=3 && user_offset~=4)
        error('User input not valid.')
    end
end

prompt3 = sprintf('\nWhich operating system are you using?\nEnter 1 or 2 to select your choice. \n 1: Windows \n 2: Mac or Linux \n');
user_OS=input(prompt3);
if(user_OS ~= 1 && user_OS ~= 2)
    error('User input not valid.')
end




%% inputs that depend on dataset
%%%%%%%%%%%%%%%%%%% GREENLAND 2019 %%%%%%%%%%%%%%%%%%%%%%%%
if(user_G_or_A == 1) 
    inputs = struct( ...
            'dataName', 'Greenland2019 Bistatic Bowtie', ...
            'startLat', 70.557821,  ... % South is negative
            'startLon', -50.091945, ... %west is negative
            'MAX_OFFSET', 2400, ...
            'SAMPLING_RATE', 15360000, ...
            'THICKNESS', 1028, ...
            'CENTER_FREQUENCY', 330e6,...
            'Chirp_Len', 8, ...  %seconds
            'gpsYN', 1, ...
            'dataType',   'short', ...       %data type from the SDR file
            'fileType',  '*.dat');       %file type of data files (should be .dat)
    if(user_OS ==1)
            inputs.dataPath = '../data/raw_data/Greenland_7_10_2019/E312bistaticBowtie0716Greenland';
    else
            inputs.dataPath = '..\data\raw_data\Greenland_7_10_2019\E312bistaticBowtie0716Greenland';
    end
end
         
     
if(user_G_or_A == 2)      
 %%%%%%%%%%%%%%%%%%% WHILLANS ICE STREAM, ANTARCTICA 2018-2019 %%%%%%%%%%%
    inputs = struct(...
            'startLat', NaN ,  ... % South is negative
            'startLon', NaN, ... %west is negative
            'MAX_OFFSET', 1300,  ...
            'SAMPLING_RATE', 15360000, ...
            'THICKNESS', 797, ...
            'CENTER_FREQUENCY', 330e6, ...
            'gpsYN', 0, ...
            'dataType',   'short', ...       %data type from the SDR file
            'fileType' , '*.dat'  );     %file type of data files (should be .dat)
    if (user_offset ==1)
            inputs.dataName = 'Whillans2019 300m Static Offset';
            if(user_OS==1)
                inputs.dataPath = '../data/raw_data/Antarctica_12_13_2018/dx0300m/slw-bistatic-dx0300-i132-f330';
            else
                inputs.dataPath = '..\data\raw_data\Antarctica_12_13_2018\dx0300m\slw-bistatic-dx0300-i132-f330';
            end
    end
    
    if (user_offset ==2)
            inputs.dataName = 'Whillans2019 600m Static Offset';
            if(user_OS==1)
                inputs.dataPath = '../data/raw_data/Antarctica_12_13_2018/dx0600m/slw-bistatic-dx0600-i132-f330';
            else
                inputs.dataPath = '..\data\raw_data\Antarctica_12_13_2018\dx0600m\slw-bistatic-dx0600-i132-f330';
            end
    end
    
    if (user_offset ==3)
            inputs.dataName = 'Whillans2019 900m Static Offset';
            if(user_OS==1)
                inputs.dataPath = '../data/raw_data/Antarctica_12_13_2018/dx0900m/slw-bistatic-dx0900-i132-f330';
            else
                inputs.dataPath = '..\data\raw_data\Antarctica_12_13_2018\dx0900m\slw-bistatic-dx0900-i132-f330';
            end
    end
    
    if (user_offset ==4)
            inputs.dataName = 'Whillans2019 1300m Static Offset';
            if(user_OS==1)
                inputs.dataPath = '../data/raw_data/Antarctica_12_13_2018/dx1300m/slw-bistatic-dx1300-i132-f330';
            else
                inputs.dataPath = '..\data\raw_data\Antarctica_12_13_2018\dx1300m\slw-bistatic-dx1300-i132-f330';
            end
    end
    
end

if ~exist(inputs.dataPath, 'dir')
    disp('File path to data does not match requested file structure. \nPlease select folder that holds your data/');
    inputs.dataPath = uigetdir;
end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%
%inputs that don't depend on which dataset
inputs.C = 299792458 ; %speed of light (m/s)
inputs.ER_ICE = 3.18; %ice permittivity
inputs.UPSAMPLE_FACTOR = 10; %upsample factor for interpolation
inputs.corr_scale_factor = 10000000; %factor for scaling to avoid rounding errors
inputs.SDR_RESISTANCE = 50; %impedance at SDR for calculating power (ohms)
inputs.user_OS=user_OS;

save('../data/processing_inputs.mat', 'inputs')