% 6/6/2019
% Nicole Bienert
% Purpose: To parse GPS data and radar data produced by run_me_GPS.sh from
% the E312 and plot the freq+time domain of the data.
% Data is formatted as ~8s of Voltage data then 100 bytes of GPS data.
% ver1: modified form of readData_ver3 for E312

%Instructions: change the variable "directory" to point towards the .dat
%file produced by the E312 SDR. If you modified the inputs to
%run_me_GPS.sh, then you may need to modify fs, dataType, and/or centerF. 

clc 
close all
clear all

directory = '/media/radioglaciology/USB31FD/test7/Thwaites2019_E312_wGPS_freq320_gain62_BW15360000_RecordLen4_SN310725C_0.dat';

fileType = '*.dat';

fs = 15.36e6; %sampling rate
T=1/fs;
dataType='short'; %data type written by SDR
centerF=320e6; %center frequency

addpath('../../myMatlabTools')
addpath('..\..\myMatlabTools')

%% Read Data

%read data
fID = fopen(directory); %open data file

data = single(fread(fID,'short'));
%100 bytes long=800 bits of GPS data 
status=fseek(fID,-100,'eof') %moves to 100 bytes from end of file
gpsData=fread(fID,'ubit8'); %read the last 100 bytes
fclose(fID);
data =data(1:length(data)-50); %The last 50 shorts are the GPS data
data = complex(data(1:2:length(data)),data(2:2:length(data))); %data is complex

%turn GPS data into a string
gpsData=char(gpsData);
gpsData=string(reshape(gpsData,[100,1])');


%% Process Data
t = (0:length(data)-1)*T;
f=linspace(centerF-fs/2,centerF+fs/2,length(data));

%perform fft
freqDomainSig = fftshift(fft(data))/length(data);
freqDomainSigMag = 20*log10(abs(freqDomainSig));

%plot freq domain
figure()
plot(f/1e6,freqDomainSigMag)
hTitle = title('Freq Domain')
hXLabel = xlabel('Freq (MHz)')
hYLabel = ylabel('20log(Quant Num) (dB)')
Aesthetics_Script

%plot time domain
figure()
plot(t,abs(data))
hTitle = title('Raw Time Domain Data')
hXLabel = xlabel('Time (s)')
hYLabel = ylabel('Quantization Number')
Aesthetics_Script