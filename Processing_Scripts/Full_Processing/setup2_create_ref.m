% 2/24/19
% Nicole Bienert
%Purpose: Read binary file for reference chirp, finds the chirp through 
%autocorrelation and exports the .mat file of the reference

clc 
close all
clear all


%% Vars
file = '/media/radioglaciology/USB31FD/test7/Thwaites2019_E312_wGPS_freq320_gain62_BW15360000_RecordLen4_SN310725C_0';
fileType = '.dat';
fs = 15.36e6;
chirpLen=1; % chirp length in seconds
sampsPerChirp=chirpLen*fs;
dataType='short';
%center freq = 320
%dur = 4s
%gain = 62

T=1/fs;
shift=1000;
%% Read Data
fID = fopen([file,fileType]); %open data file
rawDataRead = single(fread(fID,dataType));
fclose(fID);
data = complex(rawDataRead(1:2:length(rawDataRead)),rawDataRead(2:2:length(rawDataRead))); %data is complex
%use only part of the data for processing, so that the code is faster
t = (0:length(data)-1)*T;


%% Process Data
%make a rectangle to correlate against

rect=[zeros(1,shift),ones(1,sampsPerChirp),zeros(1,shift)];
rectT= [0:T:(length(rect)-1)*T];

% perform correlation with rectangle
[R,lags] = xcorr(data,rect);

% take only 1 side
R = R(lags>=0);
lags = lags(lags>=0)/fs;



%Calculate center of chirp
%correct for the zero buffer in the rectangle
Rcorrect=[zeros(shift,1);R];
%Normalize data to R and find where R is max in the chirp
%%
[middleChirpInd,val]=peakseek(abs(Rcorrect),200000,max(abs(Rcorrect))/2);
% max(real(Rcorrect(1:length(data))).*abs(data)*max(real(Rcorrect))/max(abs(data)));

%I don't care about phase, so I'm only using the real part of R
figure()
hold on
plot(t,abs(data)*max(abs(R))/max(abs(data)))
hold on
plot(lags+shift*T,abs(R),'r')
hold on
scatter(t(middleChirpInd),abs(Rcorrect(middleChirpInd)),'o')
text(t(middleChirpInd)',double(abs(Rcorrect(middleChirpInd)))',cellstr(num2str([1:length(middleChirpInd)]')));
hTitle = title('Data and xCorr with Rectangle')
hYlabel= ylabel('Magnitude')
hXlabel=xlabel('Delay (sec)')
legend('Data','xCorr with Rect')
Aesthetics_Script

userInput=input('Which matched filter peak best aligns with the center of your chirp? \nSelect your choice by entering the corresponding number.\n')
requestedChirp=middleChirpInd(userInput);
%grab chirp from time domain
x1 = requestedChirp - sampsPerChirp/2; 
x2 = requestedChirp + sampsPerChirp/2;
% take subset of data
ref_t = t(x1:x2);
ref_data = data(x1:x2);

figure()
plot(ref_t,abs(ref_data))
hTitle = title('Reference Chirp Saved')
hYlabel= ylabel('Magnitude')
hXlabel=xlabel('Delay (sec)')
Aesthetics_Script

ref_chirp=ref_data;
save('ref_chirp_BW15-36MHz_fc320_chirpLen1s.mat','ref_chirp')