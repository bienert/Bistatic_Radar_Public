%author: Nicole Bienert, Rohan Sanda
%Purpose: Apply fine scale alignment to match filtered data, coherently
%         sum and average phase-aligned chirps in same Fresnel bin.
%date: 8/31/2021


clc
close all
clear all

%% User Input
% ------------------ CONTROLS ------------------------
%for windows OS
addpath('../functions')
addpath('../functions/natsortfiles')
addpath('../data')
addpath('../data/processed_data/example_data_for_phase_2')
%for unix OS
addpath('..\functions')
addpath('..\functions\natsortfiles')
addpath('..\data')
addpath('..\data\processed_data\example_data_for_phase_2')



 filename= '../data/processed_data/';   %where to save all plots
 
 %trouble-shooting controls
 mod_freq=    10;                               %interval to save bin plots
 save_figs=   false;                             %save plots?
 addEnvelope= false;                            %applies enveloping
% ------------------ END CONTROLS ----------------------

%Load data from phase 1
load allMatchFiltChirps.mat                 
load antennaSeparation 
load processing_inputs.mat                    %input parameters from experiment

allMatchFiltChirps= allMatchFiltChirps/10^10; %make values more manageable

%------------------- CONSTANTS -------------------------                            
SIGNAL_FREQUENCY=  inputs.CENTER_FREQUENCY;      %hz
ICE_THICKNESS=     inputs.THICKNESS;             %m
ER_ICE=            inputs.ER_ICE;                %dimensionless          
C=                 inputs.C;                     %m/s
SAMPLING_RATE=     inputs.SAMPLING_RATE;         %hz 
UPSAMPLE_FACTOR=   inputs.UPSAMPLE_FACTOR;           
correlation_scale_factor= inputs.corr_scale_factor;     %scales up normalized
                                                     %matched filter output
if (inputs.user_OS ==1)
    output_folder='../data/processed_data/';
else
    output_folder='..\data\processed_data\';
end
%----------------- END CONSTANTS -----------------------

%%
%Derived constants
V=  C/sqrt(ER_ICE);                               %m/s (wave speed in ice)
num_fft_samps= ... 
    2^nextpow2(length(allMatchFiltChirps(1, :))); %number of FFT samps
upsampling_rate= SAMPLING_RATE*UPSAMPLE_FACTOR;   %hz (upsampling rate)                                        
delayDirPath= max(antennaSeparation)/C;           %s (predicted direct path)
delayBed= ...
    (2*sqrt((max(antennaSeparation)/2)^2 + ...
    (ICE_THICKNESS)^2))/V;                        %s (predicted bed path)
delay= delayBed-delayDirPath;                     %s  
startBuf= round(delay*SAMPLING_RATE * ...
           UPSAMPLE_FACTOR/4);                    %windowing frame start  
stopBuf= round(delay*SAMPLING_RATE * ...
          UPSAMPLE_FACTOR*1.5);                   %windowing frame stop                       

%Empty vectors
 numKeptChirps=        []; %vector containing the number of kept chirps per bin
 discardedCounter=     []; %vector containing the number of discarded chirps per bin
 keptCounter=          []; %vector containing the number of kept chirps per bin
 numChirpsInBin=       []; %vector containing the number of total chirps per bin
 correlationCounter=          []; %vector containing the average correlation of chirps per bin
 correlationCounterSD=        []; %vector containing the std(correlation) of chirps per bin
 meanPhaseDifference=  []; %vector containing the average hase difference of chirps per bin
 stdPhaseDifference=   []; %vector containing the std(correlation) of chirps per bin
 
 disp('STARTING PHASE 2')
 %Apply parabolic-time-shifting to matched filtered chirps
 allMatchFiltChirps= parabolicTimeShift(SAMPLING_RATE, allMatchFiltChirps);
 disp('Completed Parabolic Time Shifting')
 
%% fine scale align and coherently sum the matched filtered chirps
for ind_offset = 1:length(antennaSeparation)
    disp(['INDEX: ' num2str(ind_offset)])    
    finalInd = 1;                                                                     
    offset = antennaSeparation(ind_offset);                                           
    
    %This method ensures that the sinc functions do not seperate by more
    %than the range resolution. maxSepInBin refers to the maximum 
    %antenna separation centered around the current offset. These bounds
    %are used to find chirps in the same Fernel zone. 
    maxSepInBin = sqrt(4 * (V / (2 * SAMPLING_RATE) + ...
                  sqrt(offset^2 / 4 + ICE_THICKNESS^2))^2 - ...
                  4 * ICE_THICKNESS^2);       
    minSepInBin = sqrt(4 * (sqrt(offset^2 / 4 + ...
                  ICE_THICKNESS^2) - V / (2 * SAMPLING_RATE))^2 - ...
                  4 * ICE_THICKNESS^2);       
    chirpsInBin = find(antennaSeparation > minSepInBin & ...
                  antennaSeparation < maxSepInBin);   
    
    %Delay predictions (s)
    delayDirPath = offset / C;                             
    delayBed = (2 * sqrt((offset / 2)^2 + (ICE_THICKNESS)^2)) / V;    
    delay = delayBed - delayDirPath;                        

    %match filter against a fake direct path and basal reflection to align
    %center chirp
    t = [0:(1 / upsampling_rate):(delay * 2)];            %(s)
    
    %construct center chirp to align to
    centChirp = 4 * gauss(t, 2.1701e-08, delay / 4);      %fake direct path
    centChirp = centChirp + ...
                gauss(t, 2.1701e-08, delay / 4 + delay);  %add bed reflection
            
    %upsample in frequency domain
    centChirp_fft_up = fft(centChirp, num_fft_samps * UPSAMPLE_FACTOR);     
    centChirp_fft_up_cnj = conj(centChirp_fft_up);      

    %fft and upsample real chirp
    centChirpInd = ind_offset;
    chirp_fft = fft(abs(allMatchFiltChirps(centChirpInd, :)), num_fft_samps);   
    
    %zero pad to upsample with sinc interpolation (LPF)
    chirp_fft_up = [ ...
        zeros(1, (num_fft_samps * (UPSAMPLE_FACTOR-1)) / 2), ...
        fftshift(chirp_fft), ...
        zeros(1, (num_fft_samps * (UPSAMPLE_FACTOR-1)) / 2) ...
        ];
    
    %normalize and inverse fourier transform back to time domain
    chirp_fft_up = ifftshift(chirp_fft_up) * UPSAMPLE_FACTOR; 
    %compute upsampled chirp
    chirp_up = ifft(chirp_fft_up);
    
    
    %match filter with center chirp to find direct path peak
    R = ifft(chirp_fft_up .* centChirp_fft_up_cnj); 
    [pks pkInd] = max(abs(R));                  
    indDirect = pkInd + round(delay / 4 * ...
                (SAMPLING_RATE * UPSAMPLE_FACTOR));     

    
    %perform further windowing of chirp 
    cropStart = indDirect - startBuf;                                
    cropStop = indDirect + stopBuf;                                     
    chirp_aligned_up(centChirpInd, :) = chirp_up(cropStart:cropStop);  

    %align real chirp using direct path peak index
    centChirp=chirp_up(pkInd:pkInd + length(centChirp));

    centChirp_fft_up = UPSAMPLE_FACTOR * ...
                       fft(centChirp, num_fft_samps * UPSAMPLE_FACTOR);          
    centChirp_fft_up_cnj = conj(centChirp_fft_up);            
    
    centChirp_up = ifft((centChirp_fft_up));
    
    binPhaseMean = [];
	
	disp('End center chirp align')

	disp(['BIN SIZE' num2str(length(chirpsInBin))])
    %loop through each chirp in bin (fresnel zone) and align to center 
    %chirp
    for ind = 1:length(chirpsInBin) 
        chirpInBinOffset = antennaSeparation(chirpsInBin(ind));
           %Plot venter chirp vs. bin chirp being aligned 
           gcf1=figure(1)
           subplot(2, 1, 1)
           plot(abs(allMatchFiltChirps(centChirpInd, :)))
           hold on
           plot(abs(allMatchFiltChirps(chirpsInBin(ind), :)))
           hold off
           title('Center chirp vs current chirp being processed');
       if (mod(ind_offset, mod_freq) == 0 && save_figs) 
           saveas(gcf1, fullfile(filename, ...
                  ['MainFigure1_' num2str(ind_offset) '.png']));
       end
 
        %fft and upsample current chirp
        chirp_fft = fft(abs(allMatchFiltChirps(chirpsInBin(ind), :)), ...
                    num_fft_samps); 
                
        %zero pad to upsample with sinc interpolation (LPF)
        chirp_fft_up = [ ...
                       zeros(1, (num_fft_samps * ...
                       (UPSAMPLE_FACTOR - 1)) / 2), ...
                       fftshift(chirp_fft), ...
                       zeros(1, (num_fft_samps * ...
                       (UPSAMPLE_FACTOR - 1)) / 2) ...
                       ];
                   
        %normalize and ifft
        chirp_fft_up = ifftshift(chirp_fft_up) * UPSAMPLE_FACTOR;   
        
        %compute upsampled chirp
        bin_chirp_up = ifft(chirp_fft_up);
        
        %add envelope if requested
        if (addEnvelope)
            bin_chirp_up_copy = bin_chirp_up;
            [high, lo] = envelope(abs(bin_chirp_up_copy), 40, 'peak');
            envelopeChirp = abs(high);
            sampChirp_fft = ifftshift(fftshift(fft(envelopeChirp)));
            sampChirp_up = ifft(sampChirp_fft);
            
        else 
            sampChirp_fft = chirp_fft_up;
            sampChirp_up = bin_chirp_up;
        end  
        
        %match filter with center chirp and calculate correlation
        R = ifft(sampChirp_fft.*centChirp_fft_up_cnj);
        [correlation_pre(ind) indDirect] = max(abs(R) / ...
                                    (trapz(abs(sampChirp_fft)) * ...
                                    trapz(abs(centChirp_fft_up))) * ...
                                    correlation_scale_factor);   
            
        %calculate phase shift for direct path
        phaseShift = angle(R(indDirect));            
        
        %calculate index of direct path and perform further cropping
        indDirect = indDirect + round(delay / 4 * (SAMPLING_RATE * 10));  
        cropStart = indDirect - startBuf;                                       
        cropStop = indDirect + stopBuf;                                    
       
        %Ensure window is feasible, then time align and phase align chirp
        if cropStart > 1 && cropStop < length(sampChirp_up)                         
            chirp_aligned_up_pre = bin_chirp_up(cropStart:cropStop) * ...
                                   exp(-1i * phaseShift); 

            %prevent duplicates 
            if ind == 1
                chirp_aligned_up(finalInd, :) = chirp_aligned_up_pre;
                correlation(finalInd) = correlation_pre(ind);
                finalInd = finalInd + 1;
            else
                if max(abs(abs(chirp_aligned_up_pre) - ...
                   abs(chirp_aligned_up(finalInd - 1, :)))) > std(abs(chirp_aligned_up(finalInd - 1, :)))/1000;
               
                    chirp_aligned_up(finalInd, :) = chirp_aligned_up_pre;
                    correlation(finalInd) = correlation_pre(ind);

                    subplot(2, 1, 2)
                    plot(abs(chirp_aligned_up(centChirpInd, :)))
                    hold on
                    plot(abs(chirp_aligned_up(finalInd, :)))
                    hold off 

                    finalInd = finalInd + 1;
                end
            end
            pause(0.00001)    
         end
        
        binPhaseMean(end + 1) = ((2 * pi * SIGNAL_FREQUENCY) / C) * ...
                                (-offset + 2 * sqrt(ICE_THICKNESS^2 + ...
                                (offset / 2)^2) + chirpInBinOffset - ...
                                2 * sqrt(ICE_THICKNESS^2 + ...
                                (chirpInBinOffset / 2)^2));  
    end
    
	disp('Completed bin alignment')	

    numChirpsInBin(end + 1) = length(chirpsInBin);
    correlationCounter(end + 1) = mean(correlation);
    correlationCounterSD(end + 1) = std(correlation);
    meanPhaseDifference(end + 1) = mean(binPhaseMean);
    stdPhaseDifference(end + 1) = std(binPhaseMean);

    correlation(1) = mean(correlation);
    
    gcf3=figure(3)
    y = correlation;   
    bar(y);
    holder = [mean(correlation)];
    holder(end + 1) = std(correlation);
    title('correlation of bin');
    xlabel('Amplitude of correlation (V*sqrt(s))')
    ylabel('Each Chirp (index) in bin')
    subtitle(num2str(length(chirpsInBin)));
    subtitle(num2str(holder));
    if (mod(ind_offset, mod_freq) == 0 && save_figs) 
        saveas(gcf3,fullfile(filename, ...
            ['Bar_graph_of_correlations_' num2str(ind_offset) '.fig']));
    end
    
    %ADJUST CUTOFF THRESHOLDS HERE!!
    outliers= find(correlation<mean(correlation)-2*std(correlation) | correlation > mean(correlation)+2*std(correlation));

    
    gcf4=figure(4)
    histogram(y, 7)
    subtitle(num2str(ind_offset));
    xlabel('correlation')
    ylabel('Number of bin chirps')
    if (mod(ind_offset, mod_freq) == 0 && save_figs) 
        saveas(gcf4,fullfile(filename, ...
            ['Histogram_of_correlations_' num2str(ind_offset) '.fig']));
    end
    
    gcf5=figure(5)
    for ind = 1:length(outliers)
        plot(abs(chirp_aligned_up(outliers(ind),:)))
        hold on
    end
    title(['Discarded Chirps: ' num2str(outliers)]) 
    discardedCounter(end + 1) = length(outliers);
    if (mod(ind_offset, mod_freq) == 0 && save_figs) 
        saveas(gcf5,fullfile(filename, ...
            ['Discarded_' num2str(ind_offset) '.png']));
    end
    hold off
	disp(['DISCARDED=' num2str(length(outliers))]) 
    kept = find(~(correlation <= mean(correlation) - 2 * std(correlation) | ...
           correlation >= mean(correlation) + 2 * std(correlation)));
    
    numKeptChirps(end + 1) = length(kept);
    disp(['KEPT= ' num2str(length(kept))])
    coherentSum(ind_offset, :) = sum(chirp_aligned_up(kept, :), 1);
    coherentAverage(ind_offset, :) = coherentSum(ind_offset, :) ./ ...
                                     length(kept);
      
    gcf6=figure(6)
    for ind = 1:length(kept)
        title('Kept Chirps')
        plot(abs(chirp_aligned_up(kept(ind),:)))
        hold on
    end
    keptCounter(end + 1) = length(outliers);
    title(['Kept Chirps ' num2str(length(kept))])
    hold off
    if (mod(ind_offset, mod_freq) == 0 && save_figs)
        saveas(gcf6,fullfile(filename, ...
            ['Kept_' num2str(ind_offset) '.png']));
    end
   
    clear kept
    clear chirp_aligned_up
    clear chirp_aligned_up_pre
    clear centChirp
    clear centChirp_fft_up
    clear centChirp_fft_up_cnj
    clear chirp_fft
    clear chirp_fft_up
    clear chirp_up
    clear chirpsInBin
    clear correlation 
    clear correlation_pre
    
end

%%  ------ Get Summed Kept Chirps ------


if (addEnvelope)
    save([output_folder,'CoherentlySummedEnvelope.mat'], 'coherentSum');
    save([output_folder,'CoherentlyAveragedEnvelope.mat'], 'coherentAverage');
    save([output_folder,'NumberofKeptChirpsEnvelope.mat'], 'numKeptChirps');
else
    save([output_folder,'CoherentlySummed.mat'], 'coherentSum');
    save([output_folder,'CoherentlyAveraged.mat'], 'coherentAverage');
    save([output_folder,'NumberofKeptChirps.mat'], 'numKeptChirps');
end

% -------------------------------------
len = length(antennaSeparation);

gcf7=figure(7)
x = 1:1:len;
y1 = keptCounter;
y2 = discardedCounter;
y3 = numChirpsInBin;
plot(x,y1,'b',x,y2,'b--o', x,y3, 'k');
title('Number of kept and discarded chirps vs. antenna separation index');
xlabel('antennaSeparation Index');
ylabel('number of chirps');
legend('Kept Chirps', 'Discarded Chirps', 'Total Chirps in Bin');
saveas(gcf7, fullfile(filename, ['Kept_Counter' '.png']));

gcf8=figure(8)
x = 1:1:len;
y1 = correlationCounter;
y2 = correlationCounterSD;
plot(x,y1,'r',x,y2,'b');
title('correlation Variation vs. Antenna Sep');
subtitle('Red is mean correlation, blue is std correlation');
xlabel('antennaSeparation Index');
saveas(gcf8, fullfile(filename, ['correlation_' '.png']));


gcf9=figure(9) 
x = 1:1:len;
y1 = meanPhaseDifference;
y2 = stdPhaseDifference;
plot(x, y1, 'b', x, y2, 'r');
title('Phase Difference (rad) vs. antenna separation index');
subtitle('Red is mean difference per bin, blue is std difference for each bin');
xlabel('antennaSeparation Index');
ylabel('phase difference');
saveas(gcf9, fullfile(filename, ['Phase_' '.png']));


gcf10=figure(10)
[X,Y] = meshgrid(antennaSeparation(1:length(coherentSum(:, 1))), ...
        0:1 / SAMPLING_RATE / UPSAMPLE_FACTOR * ...
        1e6:(length(coherentSum(1, :)) - 1) * 1 / ...
        SAMPLING_RATE / UPSAMPLE_FACTOR * 1e6);
    
s = surf(X, Y, 10 * log10(abs(coherentSum./max(coherentSum')')'))
s.EdgeColor = 'none';
view([0 90]);
hTitle=title(['Coherent Sum']);
hXlabel = xlabel('Distance (m)');
hYlabel =ylabel('Travel Time (\mu s)');
ax = gca;
ax.YDir = 'reverse'
Aesthetics_Script;
colorbar
% caxis([-25 0])
cmocean('thermal')
xlim([min(antennaSeparation) max(antennaSeparation)])
ylim([0 (length(coherentSum(1, :)) - 1) * 1 / SAMPLING_RATE / ... 
      UPSAMPLE_FACTOR * 1e6])
saveas(gcf9, fullfile(filename, ['Radargram_' '.png']));

%% Functions
function Y = gauss( X , S , M )
    Y = exp(-(X - M).^2 ./ S.^2) ./ (sqrt(2 * pi) .* S);
end
    
    
    
function allMatchFiltChirps_shifted = parabolicTimeShift(SAMPLING_RATE, ...
                               allMatchFiltChirps)
                           
    dims = size(allMatchFiltChirps);
    allMatchFiltChirps_shifted = zeros(dims(1), dims(2));
        
    for chirp = 1:dims(1)
        currChirp = abs(allMatchFiltChirps(chirp, :));
 
        origChirp = allMatchFiltChirps(chirp,:);
            
        [value n] = max(abs(allMatchFiltChirps(chirp, :)));
        y0 = (currChirp(n - 1));
        y1 = (value);
        y2 = (currChirp(n + 1));
           
        y=[y0; y1; y2];
        x=[(n-1)^2 n-1 1; n^2 n 1; (n+1)^2 n+1 1];
        
        sol=inv(x)*y;
        a=sol(1);
        b=sol(2);
        c=sol(3);
             
        delta_t = -1 * ((y2 - y0) / (2 * y0 - 4 * y1 + 2 * y2));

        f = linspace(-SAMPLING_RATE/2, SAMPLING_RATE/2, length(origChirp));
            
        currChirpMod = ifft(fftshift(fft(origChirp)) .* exp(1i * (2 * pi * f) * ...
                       delta_t / SAMPLING_RATE));
                              
        t1=linspace(n-2,n+2,101);
        ynow=a*t1.^2+b.*t1+c;
        
        pkPt=n+delta_t;
        xPts=[n-20:n+20];
  
%         UNCOMMENT TO SHOW PARABOLIC FITTING!        
%         figure(1)
%         plot(t1,ynow,'--')
%         hold on
%         scatter(xPts,abs(currChirp(xPts)),'ro')
%         scatter(xPts+delta_t,abs(currChirpMod(xPts)),'r*')
%         scatter(pkPt,a*pkPt^2+b*pkPt+c)
%         xlim([xPts(1) xPts(end)])
%         legend('parabola','unshifted','shifted','pk')
%         hold off
        
        allMatchFiltChirps_shifted(chirp, :) = currChirpMod;
    end
end




