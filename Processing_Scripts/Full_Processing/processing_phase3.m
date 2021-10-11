% author: Nicole Bienert, Rohan Sanda
% date: 6/28/21
% Purpose: Convert Phase 2 output to power, then find valid basal echo
%          peaks and estimate ice-depths at each antenna separation.

clc
close all
clear all

%% User input
%----------- CONTROLS ----------------
%tell matlab where to find functions ad data
addpath('../functions')
addpath('../functions/natsortfiles')
addpath('../data')
addpath('../data/processed_data/example_data_for_phase_3')
%for unix OS
addpath('..\functions')
addpath('..\functions\natsortfiles')
addpath('..\data')
addpath('..\data\processed_data\example_data_for_phase_3')

load CoherentlyAveraged 
load ref_chirp
load antennaSeparation
load processing_inputs.mat

%controls whether to apply error_checking (TRUE) or not (FALSE)
error_check = true;
%maximum allowed time difference from previous chirp (used in error check)
TIME_THRESHOLD = 2;%2;                     %us
%minimum allowed power intensity for basal echo peak (used in error check)
NOISE_THRESHOLD = 11;%11.76;                 %dB

%control for showing figures
display = true;
%-------- END CONTROLS ----------------
%---------- CONSTANTS -----------------
%Constants used for delay estimation
ICE_THICKNESS = inputs.THICKNESS;                        %m
SDR_RESISTANCE = inputs.SDR_RESISTANCE;                  %ohms
ER_ICE = inputs.ER_ICE;                           
C = inputs.C;                                            %m/s
SAMPLING_RATE = inputs.SAMPLING_RATE;                    %hz 
UPSAMPLING_FACTOR = inputs.UPSAMPLE_FACTOR;

%---------- END CONSTANTS -------------

if (inputs.user_OS ==1)
    output_folder='../data/processed_data/';
else
    output_folder='..\data\processed_data\';
end

%% Derived constants and variables
%convert from dB to linear units
NOISE_THRESHOLD=10^(NOISE_THRESHOLD/10);

%EM wave speed
V = C / sqrt(ER_ICE);                    %m/s

%chirps have already been aligned, summed, and averaged
newSummedChirps = coherentAverage;

%create empty vector to hold power
powerChirps = [];

%Increments how many chirps have not been kept
notkeptChirps = 0;

%Vector for holding snr values for valid chirp's each basal echo peak
snrs = [];

%x-pos holding vector
xPosVec = [];

%depth-estimates holding vector
ice_depths = [];

%antenna separations of kept chirps
antennaSepKept= [];

%the basal echo index of the previous chirp 
prev_ind = NaN;

%Calculate the time duration of the reference chirp
ref_chirp_dur = length(ref_chirp) / SAMPLING_RATE;

%Correct sampling rate to account for upsampling in phase 2
SAMPLING_RATE = SAMPLING_RATE * UPSAMPLING_FACTOR;

%convert time boundary to samples
TIME_THRESHOLD = TIME_THRESHOLD * 1e-6 * SAMPLING_RATE;

%% Processing

%find power by squaring the magnitude of the summed Chirps and dividing by
%the reference chirp duration.
newSummedChirps = abs(newSummedChirps.^2) / (ref_chirp_dur * SDR_RESISTANCE);

%find basal peak and add it to powerChirps vector
for ind = 1:size(newSummedChirps,1)
    currChirp = newSummedChirps(ind, :);
                                    
    %Compute the basal echo index and power using quadratic LS estimation
    [basal_echo_power basal_echo_ind depth prev_ind snr error] = ...
        getPower(error_check, currChirp, ind, prev_ind, TIME_THRESHOLD, ...
                 SAMPLING_RATE, NOISE_THRESHOLD, V, C, ER_ICE, ...
                 ICE_THICKNESS, antennaSeparation(ind), display);
     
    if (error) 
        notkeptChirps=notkeptChirps+1;         
    else
        snrs(end+1) = snr;
        ice_depths(end+1) = depth;
        powerChirps(end+1) = basal_echo_power;
        antennaSepKept(end+1)=antennaSeparation(ind);
        xPosVec(end+1) = basal_echo_ind;
    
    end
end
        
%scatter plot of power of basal peak vs. along track distance.
gcf2 = figure(2);
    %hTitle=title('Basal Reflection Peak Power vs. Along Track Distance')
    subplot(3,1,1)
    scatter(antennaSepKept ./ 2, powerChirps)
    hXlabel = xlabel('Along Track Distance (m)');
    hYlabel = ylabel('Normalized Basal Echo Power (W)');  
    Aesthetics_Script

%scatter plot of estimated ice depth vs. along track distance
    subplot(3,1,2);
    %hTitle=title('Ice thickness vs. along track distance')
    scatter(antennaSepKept ./ 2, ice_depths)
    hXlabel = xlabel('Along Track Distance (m)');
    hYlabel = ylabel('Estimated Ice Depth (m)');  
    set(gca, 'YDir','reverse')
    Aesthetics_Script

 %Plot of basal echo power vs. path length with one way average attenuation
 %coefficient    
    subplot(3,1,3);
    pathLen=2*sqrt((antennaSepKept ./ 2).^2 + ice_depths.^2);
    amp=10*log10(powerChirps ./ (C ./ (4 * pi * SAMPLING_RATE ... 
                                            * pathLen / (sqrt(ER_ICE)))).^2);
    p=polyfit(pathLen,amp,1);
    attenAvg=p(1)*1000; %Two-way avg atten
    scatter(pathLen,amp)
    hold on
    plot(pathLen,polyval(p,pathLen))
    hold off
    xl=xlim;
    yl=ylim;
    text(mean(xl), mean(yl) + (yl(2) - yl(1))/4, ...
        ['One Way Avg Atten = ', num2str(attenAvg)], ...
        'FontName', 'Times New Roman', 'FontSize', 11);
    hXlabel = xlabel('Path Length (m)');
    hYlabel = ylabel('Normalized Basal Echo Power (dB)'); 
    Aesthetics_Script
    
disp('Done!!')

%save basal echo power and estimated ice depths
save([output_folder,'basalEchoPowers.mat'], 'powerChirps');
save([output_folder,'ice_depths.mat'], 'ice_depths');


%% Function Definitions
%This function checks whether a given chirp has a feasible basal echo peak,
%or whether it requires further processing. An output boolean of TRUE 
%indicates that the peak is an error, and not feasible. 
function [error prev] = checkBasalPeak(error_check, ind, prev, pos, ...
                                        TIME_THRESHOLD, delayBed, SAMPLING_RATE,...
                                        NOISE_THRESHOLD, snr);
    if (~error_check)
        error = false;
        prev = NaN;
        return;
    end
    
    if ind == 1 || prev == NaN
        delaySamps = delayBed * SAMPLING_RATE;
        if abs(pos - delaySamps) > TIME_THRESHOLD || snr < NOISE_THRESHOLD
            error = true;
        else 
            error = false;
            prev = pos;
        end  
    else
        if abs(pos - prev) > TIME_THRESHOLD || snr < NOISE_THRESHOLD
            error = true;
        else
            error = false;
            prev = pos;
        end
    end
end

%This function first finds the peak amplitude and index of the basal echo
%peak. It then fits a parabola to three points of the echo peak using 
%the Quadratic LS Peak Estimation method to recover the true amplitude 
%lost during discrete interpolation. 
function [power peak_ind depth prev_ind snr error] = ...    
                                  getPower(error, chirp, ind, prev,...
                                           TIME_THRESHOLD, SAMPLING_RATE,...
                                           NOISE_THRESHOLD, V, C, ER_ICE,...
                                           ICE_THICKNESS, distance, disp)
                                      
    currChirp = chirp;
    origChirp = chirp;

    %Find the peak amplitude and index
    basalHalfInds = round(length(currChirp) / 2):length(currChirp);
    basalHalf = currChirp(basalHalfInds);

    [value index] = max(basalHalf);
    n = basalHalfInds(index);

    %parabolic fitting
    y0 = (currChirp(n - 1));
    y1 = (value);
    y2 = (currChirp(n + 1));

    %Quadratic LS Peak Estimation
    delta_t = -1 * ((y2 - y0) / (2 * y0 - 4 * y1 + 2 * y2));
    correctIndex = n - delta_t;

    y = [ y0;    y1; y2];
    x = [ n - 1; n;  n + 1 ];
    B = [ x.^2   x   ones(3, 1) ]; 

    coeffs = y\B;

    a = coeffs(1);
    b = coeffs(2);
    c = coeffs(3);

    %resample waveform to new peak
    f = linspace(-SAMPLING_RATE/2, SAMPLING_RATE/2, length(origChirp));  
    currChirpMod = ifft(fftshift(fft(origChirp)) .* exp(1i * (2 * pi * f) * ...
                   delta_t / SAMPLING_RATE));  

    currChirpMod = abs(currChirpMod);
    basal_power = currChirpMod(n);
    
    %estimate snr of basal peak
    basalHalfIndsMod = round(length(currChirpMod) / 2):length(currChirpMod);
    basalHalfMod = currChirpMod(basalHalfIndsMod);
    noise=[basalHalfMod(1:n-200-round(length(currChirpMod) / 2)), ...
           basalHalfMod(n+200-round(length(currChirpMod) / 2):end)];

    snr = basal_power / std(noise);
    
    %estimate the bed echo peak delay for error checking
    delayBed = (2 * sqrt((distance / 2)^2 + (ICE_THICKNESS)^2)) / V;
    
    [error prev] = checkBasalPeak(error, ind, prev, n, ...
        TIME_THRESHOLD, delayBed, SAMPLING_RATE, NOISE_THRESHOLD, snr);
      
    if (error)
        power = NaN;
        peak_ind = NaN;
        depth = NaN;
        prev_ind = prev;
        snr = snr;
        error = error;
    else
        %Depth Estimation
        tb = correctIndex / SAMPLING_RATE;
        [direct_peak direct_index] = max(currChirp);
        td = direct_index / SAMPLING_RATE;

        t0 = td - distance / C;

        h = sqrt((C * (tb - t0) / (2 * sqrt(ER_ICE)))^2 - (distance)^2 / 4);
        depth = h;

        if (disp)
            %Plot the chirp with the basal peak highlighted.
            figure(1);               
            plot(currChirpMod,'b')
            hold on
            plot(correctIndex,basal_power,'o','MarkerSize', 10);
            hold off
            yl = ylim;
            ylim(yl / 100)
            pause(0.1)

         end
         power = basal_power;
         peak_ind = correctIndex;
         prev_ind = prev;
         snr = snr;
    end
end
