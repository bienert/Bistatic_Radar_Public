% author: Nicole Bienert, Rohan Sanda
% date: 8/31/21
% Purpose: Matched filter received signal with reference chirp and apply 
%          coarse time alignment to chirps. 

%Important Notes:
%   The raw time domain data should always have a quantization number <
%   6,000. Any higher indicates that the maximum input voltage to the SDR
%   has been exceeded and the SDR will get damaged. 

clc 
close all
clear all

%for windows OS
addpath('../functions')
addpath('../functions/natsortfiles')
addpath('../data')
%for unix OS
addpath('..\functions')
addpath('..\functions\natsortfiles')
addpath('..\data')


load processing_inputs.mat  %input parameters 
load ref_chirp.mat          %load pre-recorded chirp

%% User Input
%------------- CONTROLS ------------------
display=     2;            %can be set to 0, 1, 2, or 3 
                               %adjusts what plots are shown 
save_data=   false;         %save all data files?
save_plots=  false;         %save all plots?

%paths to save match filter and raw time domain outputs (in that order).
%only applicable if display > 1
filenames= ['../phase1_output/match', '../phase1_output/raw']; 
                                                                

%------------- END CONTROLS --------------
%------------- INPUT VARIABLES ----------- 
startLat= inputs.startLat;                  %starting latitude of transmit 
                                                %antenna (use - for S)
startLon= inputs.startLon;                  %starting longitude of  
                                                %antenna(use - for W)
SAMPLING_RATE= inputs.SAMPLING_RATE;        %sampling rate
THICKNESS=     inputs.THICKNESS;            %ice thickness (m) 
MAX_OFFSET=    inputs.MAX_OFFSET;           %max separation (m)
C=             inputs.C;                    %speed of light (m/s)

dataFolder=    inputs.dataPath;             %data folder file path
gpsYN=         inputs.gpsYN;                %is GPS data embedded? (1) yes (0) no
dataType=      inputs.dataType;             %data type from the SDR file
fileType=      inputs.fileType;             %file type of data files (should be .dat)
myTitle=       inputs.dataName;             %Figure title

if (inputs.user_OS ==1)
    output_folder='../data/processed_data/';
else
    output_folder='..\data\processed_data\';
end
%--------- END INPUT VARIABLES ------------

%% Internal variables
%Derived Constants
T=                1/SAMPLING_RATE;              %sample period
sampsPerChirp=    1/200e6*...                   %number of expected samples
                  (SAMPLING_RATE)^2;                %in a chirp                                                
anticipatedDelay= 2*sqrt(THICKNESS^2+MAX_OFFSET^2) ...
                  /C-MAX_OFFSET/C;              %estimate of delay between
                                                     %bed and direct path
windowSize=       floor(anticipatedDelay * ...  %Display window for viewing
                  4/T*8/2)*2;                        %match filtered data
fresnelZone=      sqrt(THICKNESS*C/ ...         %Fresnel zone approximation
                  (2*(320e6-SAMPLING_RATE/2)));


%% Get filenames
%Create directory of filenames sorted by date
b= subdir(fullfile(dataFolder, fileType));
S= [b(:).datenum].';                            %obtain date

disp("Starting Phase 1 Processing.")
disp("Data files being processed: ")

file_names=   {b(:).name}.';                    %get file names
sorted_names= natsortfiles(file_names);         %sort file names by number
directory=    sorted_names;                                  
disp(directory)                                 %print dates for reference 

dataFolderStr= regexprep(dataFolder,'\','\\\');

%obtain the filename of each recording in the 'data' folder
name= regexprep(regexprep(regexprep(regexprep(directory, ...
                          dataFolderStr,''), '\\',''), ...
                          regexprep(fileType,'*',''),''),'_','-');

set(0,'defaultAxesFontsize',18)

%Empty vectors
allNoChirpsFiles=   [];     %vector containing file names with no chirps
numNoChirpsFiles=   0;      
totalNumChirps=     0;
allMatchFiltChirps= [];     %complex matrix containing phase 1 output
gps_chirps=         [];     %contains gps coordinates for each chirp
lastGPSdistance=    0;      %useful for trouble-shooting
noisePWs=           [];     %average signal power minus direct path peak
noiseFlrs=          [];     %average correlation + 10*std(correlation)
times=              [];     %GPS times for each data file
antennaSeparation=  [];     %transmitter+receiver offset per chirp
velocity_vec=       [];     %average velocity for each data file
gpsPts=             [];     %GPS coordinates (lat/lon) for each data file

%% Begin processing
m = 1;
for m = m:length(directory)
    disp("File being processed: ")
    disp(char(directory{m})) 
    
    [matchFiltSub, numChirpsDetected, noChirpsFiles, ...
    chirpTimeFromEnd, gpsData, noisePW, noiseFlr] = ...
        matchedFilterAlign(m, ...
            ref_chirp, SAMPLING_RATE, char(directory{m}), ...
            windowSize, dataType, display, gpsYN, filenames); 
	    
    %check if there were chirps in the file
    if ~isempty(noChirpsFiles)
        disp('REJECTED. No chirps found.')
        numNoChirpsFiles= numNoChirpsFiles+1;
        allNoChirpsFiles{numNoChirpsFiles}= noChirpsFiles;
        if numChirpsDetected == 0
            noiseFlrs(end+1)= noiseFlr;   
            noisePWs(end+1)= NaN;          
        end
    else
        %add detected chirps to array
        allMatchFiltChirps= [allMatchFiltChirps; matchFiltSub];  
        totalNumChirps= totalNumChirps+numChirpsDetected;
        %update distances for measurements
        if numChirpsDetected > 0
	       noisePWs(end+1)= noisePW;
           noiseFlrs(end+1)= noiseFlr;
        end
    end
    
    if gpsData.valid 
	    disp('VALID GPS DATA.')
        gpsPts= [gpsPts, [gpsData.lat; gpsData.lon]];
        gpsDistance= haversine([startLat, startLon], ...
                               [gpsData.lat, gpsData.lon]);
        if m == 1
            velocity= -5.777302248359114e-08;
        else
            velocity= (gpsDistance-lastGPSdistance)/ ...
                       (gpsData.GMTtime-lastGMTtime);
        end
        velocity_vec(end+1)= velocity;
        
        if isempty(noChirpsFiles)
            antennaSeparation= [antennaSeparation, gpsDistance - ...
                                velocity.*chirpTimeFromEnd];
            gps_chirps= [gps_chirps, ... 
                        [gpsData.lat.*ones(1, numChirpsDetected); ... 
                         gpsData.lon.*ones(1, numChirpsDetected)]];
        end
        
        lastGMTtime= gpsData.GMTtime;
        times(end+1)= gpsData.GMTtime;
        lastGPSdistance= gpsDistance;

    else
        disp('INVALID GPS DATA.')
	    gpsPts= [gpsPts, [NaN; NaN]];
        velocity_vec(end+1)= NaN;
	    times(end+1)= NaN;	
	    if isempty(noChirpsFiles)
	   	     antennaSeparation= [antennaSeparation, chirpTimeFromEnd*0];
	         gps_chirps=  [gps_chirps ... 
                          [zeros(1, numChirpsDetected);  ...
                           zeros(1, numChirpsDetected)]];
	    end
   	end   
                 
    %Save important matricies in case it crashes
     if ~isempty(matchFiltSub) && save_data
         save([output_folder,'allMatchFiltChirps.mat'],'allMatchFiltChirps')
         save([output_folder,'antennaSeparation.mat'],'antennaSeparation')
         save([output_folder,'gpsPts.mat'],'gpsPts')
         save([output_folder,'noisePWs.mat'], 'noisePWs')
         save([output_folder,'noiseFlrs.mat'], 'noiseFlrs')
         save([output_folder,'times.mat'], 'times')
         save([output_folder,'velocity_vec.mat'], 'velocity_vec')
	     save([output_folder,'gps_chirps.mat'], 'gps_chirps')
     elseif ~isempty(allNoChirpsFiles) && save_data
         save([output_folder,'numNoChirpsFiles.mat'],'numNoChirpsFiles')
         save([output_folder,'allNoChirpsFiles.mat'],'allNoChirpsFiles')
     end
end

save([output_folder,'allMatchFiltChirps.mat'],'allMatchFiltChirps')
save([output_folder,'antennaSeparation.mat'],'antennaSeparation')
save([output_folder,'gpsPts.mat'],'gpsPts')
save([output_folder,'noisePWs.mat'], 'noisePWs')
save([output_folder,'noiseFlrs.mat'], 'noiseFlrs')
save([output_folder,'times.mat'], 'times')
save([output_folder,'velocity_vec.mat'], 'velocity_vec')
save([output_folder,'gps_chirps.mat'], 'gps_chirps')



if numNoChirpsFiles < length(directory)
     %instrument may have moved backwards, so re-organize data
     [antennaSeparation, sortIndex]= sort(antennaSeparation, 'ascend');
     allMatchFiltChirps= allMatchFiltChirps(sortIndex, :);
 

     
     %Path plot of instrument on map
     if(gpsYN ==1)
         
          %Radargram of matched filtered chirp
         gcf4= figure(4);
         [X, Y]= meshgrid(antennaSeparation, 0:T*1e6:(windowSize)*T*1e6);
         s= surf(X, Y, 10*log(abs(allMatchFiltChirps')));
         s.EdgeColor= 'none';
         view([0 90]);
         hTitle= title(['Match Filter of ', num2str(totalNumChirps), ...
                        ' Chirps at ', myTitle]);
         hXlabel= xlabel('Distance'); 
         hYlabel= ylabel('Travel Time (microseconds)');
         ylim([25 100])
         Aesthetics_Script;
         if save_plots
             saveas(gcf, [myTitle, '_quasiRadarGram_Distance_', ...
                    num2str(totalNumChirps), '_Chirps'], 'fig')
             saveas(gcf, [myTitle, '_quasiRadarGram_Distance_', ...
                    num2str(totalNumChirps), '_Chirps'], 'png')
         end

     
         gcf5= figure(5);
         hold on
         worldmap([min(gpsPts(1, :))-0.000001 max(gpsPts(1, :))+0.000001], ...
                  [min(gpsPts(2, :))-0.000001 max(gpsPts(2, :))+0.000001])
         geoshow(gpsPts(1, :), gpsPts(2, :), 'displaytype', 'point')
         geoshow(gpsPts(1, :), gpsPts(2, :), 'displaytype', 'line')
         geoshow(gpsPts(1, :), gpsPts(2, :), 'displaytype', 'point')
         geoshow(startLat, startLon, 'displaytype', 'point')
         for k = 1:length(gpsPts(1, :))
             textm(gpsPts(1, k), gpsPts(2, k), num2str(k))
         end
         if save_plots
             saveas(gcf5, 'gpsPts', 'fig')
             saveas(gcf5, 'gpsPts', 'png')
             saveas(gcf5, [myTitle, '_worldmap_', ...
                    num2str(totalNumChirps), '_Chirps'], 'fig')
         end
     end
 
     %plot each chirp vs chirp number
     gcf6= figure(6);
     [X, Y]= meshgrid([1:totalNumChirps], 0:T*1e6:(windowSize)*T *1e6);
     s= surf(X, Y, 10 * log(abs(allMatchFiltChirps'))) 
     s.EdgeColor= 'none';
     view([0 90]);
     hTitle= title(['Match Filter of ',num2str(totalNumChirps), ...
                    ' Chirps at ',myTitle]);
     hXlabel= xlabel('Chirp Numper ~ Distance');
     hYlabel= ylabel('Travel Time (microseconds)');
     ylim([25 100])
     Aesthetics_Script;
     if save_plots
        saveas(gcf6, [myTitle, '_quasiRadarGram_', num2str(totalNumChirps), ...
                  '_Chirps'], 'fig')
        saveas(gcf3, [myTitle, '_quasiRadarGram_', num2str(totalNumChirps), ...
                   '_Chirps'], 'png')
     end
end
