% 6/11/2019
% Nicole Bienert
% Purpose: To verify that E312 is seeing chirps and recording GPS data.
%Parses GPS and SDR data then match filters it with a pre-recorded reference chirp.

% %Instructions: change the variable "directory" to point towards the .dat
%file produced by the SDR. If you modified the inputs to
%run_me_GPS.sh then you may need to  modify fs, and dataType. This can only
% be used with the E312 becaus the format for GPS data is different for the
% x310. data is formatted as 8s of voltage data then 100 bytes of GPS data.

%Warning: This code has not been tested

clc 
close all
clear all

directory = '/media/radioglaciology/USB31FD/test2/Thwaites2020_E312_wGPS_freq320_gain62_BW15360000_RecordLen8_SN310725C_0.dat';
fileType = '*.dat';

fs = 15360000;
T=1/fs;
dataType='short'; %default data type for E312 and x310

%center freq = 330
%dur = 8s
%gain = 62


addpath('../../myMatlabTools')
addpath('../refChirps')
addpath('..\..\myMatlabTools')
addpath('..\refChirps')

%% Read Data

%This pre-recorded reference chirp only applies to E312 data collected with fs=15.36MHz
load  ref_chirp_BW15-36MHz_fc320_chirpLen1s

chirpLen=1; % chirp length in seconds
sampsPerChirp=chirpLen*fs;
T=1/fs;
sampsSepPks=20; %number of samples required to make two peaks not be the same peak 

%read data
fID = fopen(directory); %open data file

data = single(fread(fID,dataType));
%100 bytes long=800 bits of GPS data 
status=fseek(fID,-100,'eof') %moves to 100 bytes from end of file
gpsData=fread(fID,'ubit8'); %read the last 100 bytes
fclose(fID);
data =data(1:length(data)-50); %The last 50 shorts are the GPS data
data = complex(data(1:2:length(data)),data(2:2:length(data))); %data is complex

%turn GPS data into a string
gpsData=char(gpsData);
gpsData=string(reshape(gpsData,[100,1])');

% perform match filtering
R = matchedFilt(data,ref_chirp);

%Plot match filter
figure()
plot([0:T:(length(R)-1)*T],abs(R));
xlim([0, (length(R)-1)*T])
hTitle=title({'Match Filtered Data';'(Operating on Voltage Data)'})
hXlabel= xlabel('Time (seconds)')
hYlabel=ylabel('Quantization Number Squared')
Aesthetics_Script

%Plot time domain
figure()
plot([0:T:(length(data)-1)*T],abs(data).^2);
xlim([0, (length(data)-1)*T])
hTitle=title('Time Domain Data')
hXlabel= xlabel('Time (seconds)')
hYlabel=ylabel('Quantization Number Squared (\propto PW)')
Aesthetics_Script
pause(0.01)
% 
% %Plot I and Q
% figure()
% hold on
% plot([0:T:(length(data)-1)*T],imag(data));
% plot([0:T:(length(data)-1)*T],real(data));
% xlim([0, (length(data)-1)*T])
% hTitle=title('Time Domain Data')
% hXlabel= xlabel('Time (seconds)')
% hYlabel=ylabel('Quantization Number')
% hLegend = legend('imag','real')
% Aesthetics_Script
% pause(0.01)
% 
% %moving average of power in 1ms windows
% movingAveragePW=movmean(abs(data).^2,0.001*fs);
% figure()
% plot([0:T:(length(data)-1)*T],movingAveragePW);
% xlim([0, (length(data)-1)*T])
% hTitle=title('Moving Average of Time Domain Data')
% hXlabel= xlabel('Time (seconds)')
% hYlabel=ylabel('Quantization Number Squared (\propto PW)')
% Aesthetics_Script
% pause(0.01)


