% author: Nicole Bienert, Rohan Sanda
% date: 8/31/21
% Purpose: A function that performs matched filtering 


function [matchFiltSub, numChirpsDetected, noChirpsFiles, ...
          chirpTimeFromEnd, gpsData, noisePW, noiseFlr]=  ...
                matchedFilterAlign(m, ref_chirp, fs, ...
                directory, windowSize, dataType, display, gpsYN, filenames)
            
    quantScaling=  1/64191899;
    sampsPerChirp= 1/200e6*fs^2;     %always even
    T=             1/fs;
    sampsSepPks=   20;               %number of samples required to make two
                                         %peaks not be considered the same peak 

    %function that reads voltage and gps data. Works for both GPS and non-GPS
        %stamped files
    [data, gpsData]= readE312Data(directory, dataType, gpsYN);
    data= data*quantScaling;

    % perform match filtering
    [R,lags]= xcorr(data, ref_chirp);
    R= R(lags >= 0);
    
    clear lags 

    %% find the direct path for phase shift data processing
    %if bed is visible, it may be higher than the direct path. If bed is
    %visible then use only the first peak as the direct path.
    %find all direct paths and bed reflections.

    % find peaks above noise floor
    noiseFlr= mean(abs(R))+10*std(abs(R));
    threshold= max(noiseFlr, max(abs(R))/4);
    [locs, peaks] = peakseek(abs(R), sampsSepPks, threshold); 
    
    filename1= filenames(1);
    filename2= filenames(2);
    %Display returned data if requested
    if display == 2 || display == 3
        %Plot matched filter output
        gcf2=figure()
        plot([0:T:(length(R)-1)*T],abs(R));
        xlim([0, (length(R)-1)*T])
        hTitle= title({'Match Filtered Data'; num2str(m)})
        hXlabel= xlabel('Time (seconds)')
        hYlabel= ylabel('$$V \sqrt(s) $$', 'interpreter', 'latex');
        xlim([0 8])
        Aesthetics_Script
        pause(0.01)
        saveas(gcf2, fullfile(filename1, ['matchfilt_' num2str(m) '.png']))
        
        %Plot raw time domain signal (8-s chirp)
        gcf3=figure()
        plot([0:T:(length(data)-1)*T],abs(data));
        xlim([0, (length(data)-1)*T])
        hTitle = title({'Time Domain Data'; num2str(m)})
        hXlabel = xlabel('Time (seconds)')
        hYlabel = ylabel('Voltage')
        xlim([0 8])
        Aesthetics_Script
        pause(0.01)
        saveas(gcf3, fullfile(filename2, ['raw_signal_' num2str(m) '.png']))
    end
    
    %only continue if there were chirps in the file
    if isempty(locs)
        matchFiltSub= [];
        numChirpsDetected= 0;
        chirpTimeFromEnd= [];
        noisePW= [];
        noChirpsFiles= directory;
        if display == 1
            %Plot of matched filter output 
            figure()
            plot([0:T:(length(R)-1)*T],abs(R));
            xlim([0, (length(R)-1)*T])
            hTitle= title('Match Filter - No Peaks Above Threshold')
            hXlabel= xlabel('Time (microseconds)')
            hYlabel= ylabel('$$V \sqrt(s) $$', 'interpreter', 'latex')
            Aesthetics_Script
        end
    else

        %remove any peaks which have been clipped or are too close to the edge 
        %for the summation algorithm. Must use while loop since locs changes size
        n = 1;
        while (n <= length(locs))
            if ((locs(n)+sampsPerChirp/2 < (sampsPerChirp+1))) 
                locs(n)= [];
                peaks(n)= [];
            elseif ((locs(n)+sampsPerChirp/2 > (length(data)-sampsPerChirp-1)))
                locs(n)= [];
                peaks(n)= [];
            else
                n= n+1; 
            end
        end

        directPathLocs= locs;
        directPathPks= peaks;
        numChirpsDetected= length(directPathLocs);

    %%
        noise= ((R));
        shift= []; 

        for q = 1:length(directPathLocs)
            if ((directPathLocs(q)-2e6-(q-1)*4e6) > 1) && ...
                    ((directPathLocs(q)+2e6-(q-1)*4e6) < length(noise))
                        noise(directPathLocs(q)-2e6-sum(shift): ...
                            directPathLocs(q)+2e6-sum(shift))= [];
                        shift(q)= 4e6;
            elseif ((directPathLocs(q)-2e6-(q-1)*4e6) < 1)
                noise(1:directPathLocs(q)+2e6-sum(shift))= [];
                shift(q)= directPathLocs(q)+2e6-sum(shift);
            elseif ((directPathLocs(q)+2e6-(q-1)*4e6) > length(noise))
                noise(directPathLocs(q)-2e6-sum(shift):end)= [];
                shift(q)= length(noise)-directPathLocs(q)-2e6-sum(shift);
            end  
        end

        %chop ends in case chirp was clipped
        noise= noise(sampsPerChirp:end-sampsPerChirp);

        if display == 3
            %Plot noise
            gcf= figure(10)
            plot(abs(noise))
            saveas(gcf,'noise','png')
        end

        %convert noise to power
        noisePW= mean(abs(noise)).^2;
        %% zoom in on each rx chirp, coarse time alignment 
        for n = 1:numChirpsDetected
            %grab time domain data centered around chirp at 2x chirp width
           
            matchFiltSub(n,:)= R(directPathLocs(n)-windowSize/4: ...
                                    directPathLocs(n)+windowSize/4*3); 

            matchFiltSub(n,:)= matchFiltSub(n,:)* ...
                                    (exp(-1i*angle(matchFiltSub(n,floor(windowSize/4+1)))));
            %perform adjustment twice to mitigate error
            matchFiltSub(n,:)= matchFiltSub(n,:)* ...
                                    (exp(-1i*angle(matchFiltSub(n,floor(windowSize/4+1)))));
            
        end
        
        disp(size(matchFiltSub))

        noChirpsFiles= [];
        chirpTimeFromEnd= (length(R)-directPathLocs)*T;

    end %end of if statement that checks that there are chirps
end % end for function
