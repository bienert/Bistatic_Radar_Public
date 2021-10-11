# Processing Chain ReadMe
Please note that this repository accompanies a paper in review titled "Large Offset Bistatic Radar Sounding of Glaciers using Processing-Based Synchronization". Please cite this paper to comply with our license when redistributing transformed versions of this work (e.g. publishing a paper on experiments which which used some code from this repo). For example: N. Bienert et al., "Post-Processing Synchronized Bistatic Radar for Long Offset Glacial Sounding," IEEE TRANS. Geosci. Remote Sens., in review.

### Table of Contents
* [Welcome](#welcome)
* [Setup](#setup)
  * [Software](#software)
  * [Subsidiary Functions](#Subsidiary-Functions)
  * [Data Files](#Data-Files)
  * [Directory Infrastructure](#directory-infrastructure)
* [Processing Chain](#processing-chain)
  * [Terminology](#terminology)
  * [Summary](#summary)
  * [Set-Up 1](#Set-Up-1)
  * [Set-Up 2](#Set-Up-2)
  * [Phase 1](#phase-1)
  * [Phase 2](#phase-2)
  * [Phase 3](#phase-3)
* [Example](#example)
* [Note on High Performance Computing](#cluster-computing)
* [Questions](#questions)

## Welcome 
Welcome to the Processing Chain ReadMe. Below you will find a comprehensive guide to using this repository. We have provided a number of example datasets to allow you to run this repository on our field data from Greenland and Antarctica.  

## Setup
**This respository is meant to be fully executable with the example dataset after completing the steps in the Setup section.** The scripts already contain default variable settings optimized for the [example dataset](#data-files). A full description of each input variable is listed below to provide guidance in processing your own data.

### Software
All scripts were written in Matlab v2020+.

### Subsidiary Functions
Download the following functions from [MathWorks](https://www.mathworks.com/), place each one into the 'functions' subfolder shown above.
* `haversine` - Download the `haversine.m` script [here](https://www.mathworks.com/matlabcentral/fileexchange/27785-distance-calculation-using-haversine-formula)
* `peakseek` - Download the `peakseek.m` script [here](https://www.mathworks.com/matlabcentral/fileexchange/26581-peakseek)
* `subdir` - Download the `subdir.m` script [here](https://www.mathworks.com/matlabcentral/fileexchange/15859-subdir-a-recursive-file-search)
* `natsortfiles` - Download and uncompress the folder containing the `natsortfiles` function [here](https://www.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort). Place the entire folder into the `tools` folder.
* 'cmocean' - Download the 'cmocean.m' script [here](https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps)

### Data Files
While phase 2 and 3 can be run with the example dataset provided in the [processed_data](https://github.com/bienert/Bistatic_Radar/tree/main/Processing_Scripts/data/processed_data) folder, phase 1 requires a raw dataset. Download the data we collected at https://www.usap-dc.org/view/dataset/601472 and place the contents in the [raw_data] folder (see directory infrastructure below). 

### Directory Infrastructure
After cloning this repository to your local computer, ensure that the following file infrastructure exists within the `Processing Scripts` folder. 
```
Full_Processing
├──processing_phase1.m
├──processing_phase2.m
├──processing_phase3.m
├──setup1_inputs_struct.m
└──setup2_create_ref

functions
├── matchedFilterAlign.m
├── readE312Data_ver2.m
├── haversine.m
├── peakseek.m
├── subdir.m
├── cmocean.m
├── natsortfiles (folder)
└── Aesthetics_Script.m

data
├── ref_chirp.mat
├── processing_inputs.mat
├── processed_data (folder)
└── raw_data
    ├── Greenland_7_10_2019
         ├── E312bistaticBowtie0716Greenland (folder)
         ├── E312bistaticRugged0716Greenland (folder)
         └── E312bistaticNoGPSRugged0715Greenland (folder)
    └── Antarctica_12_13_2018
         ├── dx0300m
              └── slw-bistatic-dx0300-i132-f330 (folder)
         ├── dx0600m
              └── slw-bistatic-dx0300-i132-f330 (folder)
         ├── dx0900m
              └── slw-bistatic-dx0300-i132-f330 (folder)
         └── dx1300m
              └── slw-bistatic-dx0300-i132-f330 (folder)


```

## Processing Chain

### Terminology

The goal of the processing chain is to recover the **bed echo peak** using coherently summation. We use information from the **direct path peak** to time and phase-align chirps. We describe these two terms below, as they are used throughout the processing guide.

<p align="center">
  <img width="559" height="263" src="https://github.com/bienert/Bistatic_Radar/blob/main/misc/Images/terminology_fig.png"> 
</p>

The diagram above shows a simple bistatic radar configuration. The transmitter sends two chirps along two paths: the **direct path** (in blue) and the **bed path** (in red). The corresponding direct and bed path peaks are shown in the accompanying plot, which graphs the matched filtered chirp versus time. 


### Summary
We have divided our processing chain into set-up steps and three processing phases. [Set-Up 1](#Set-Up-1) must be run before any other script, because it configures variables used by phase 1-3. [Set-Up 2](#Set-Up-2) creates a reference chirp used in phase 1-3, and should only be used when processing freshly collected field data. We provide the reference chirp to accompany our data, so this step is not required when working with our datasets. [Phase 1](#phase-1) matched filters the received signal with a pre-recorded reference chirp to recover the direct path peak. [Phase 2](#phase-2) upsamples the output from Phase 1 and performs fine-scale phase alignment and coherent summation to recover the bed echo peak. Phase 2 is currently configured to read example phase 1 output data, and can be run without first running phase 1. Finally, [Phase 3](#phase-3), determines the bed echo power for each received signal. Phases 3 is configured to read in example phase 2 output data, so phase 3 can be run without first running phase 1-2.

### Set-Up 1
''Required'' To ensure that these scripts can be easily used across multiple sets of experimental data, we have streamlined all scripts' input variables into a single Matlab `struct` file. Run the `setup1_inputs_struct.m` script to select your dataset and configure the code for your operating system. Descriptions of all variables for each phase are listed in [Processing Chain](#processing-chain). The script outputs a `processing_inputs.mat` file containing the `struct` object which is automatically saved to the `data` folder. 

### Set-Up 2
This step is only required for data you collect. We provide a reference chirp for use with our datasets, but if you build the radar system and conduct an experiment following our [manual](https://github.com/bienert/Bistatic_Radar/blob/main/Bienert's%20Bistatic%20Manual.pdf), you will need to create a new reference chirp to accompany your data (as described in the manual).

### Phase 1

<p align="center">
  <img width="675" height="350" src="https://github.com/bienert/Bistatic_Radar/blob/main/Misc/Images/phase1.png"> 
</p>

#### Description
Phase 1 begins by first calling the function `matchedFilterAlign` on each 8-s time-domain file from the SDR. This function first extracts GPS-position information from the last 100 bytes of each file using the `readE312Data` function. It then matched filters each time-domain file with the reference chirp (`ref_chirp`) to recover the direct path of each transmitted chirp. As shown in the Phase 1 figure, the function clips the data file around each peak, thereby isolating each chirp, and proceeds to coarsly time align the chirps. An example of an outputted chirp is shown in the last plot in the Phase 1 figure. After iterating through all the data files, the script assigns antenna separation offsets to each individual chirp using the GPS information. All chirps are stored in the **allMatchFiltChirps** file and the corresponding antenna separations are stored in **antennaSeparation**.
#### Input
All input constants (in the **Constants** section) are automatically loaded from the `processing_inputs.mat` file to ensure the continuity of values throughout all phases of processing. See the [Variable Input Methods](#variable-input-methods) section for information on how to edit the file. All input controls (in the **Controls** section) should be entered manually.
##### Controls
* **dataFolder** [array-character] - a filepath to the `data` folder containing the raw data files
* **display** [integer] - enter an integer of either 0, 1, 2, or 3 to control output plots
  * 0 = only display final matched filter output
  * 1 = additionally, display matched filtered out of data where no correlation peaks were detected
  * 2 = additionally, display the raw time-domain signal
  * 3 = additionally, display aligned chirps (trouble-shooting mode)
* **gpsYN** [integer] - enter an integer 0 or 1 to indicate gps type
  * 0 = FALSE, the gps/transmitter is stationary
  * 1 = TRUE, the gps/transmitter is moving along a transect 
* **dataType** [array-string] - the data type of the SDR file (ex. 'short')
* **fileType** [array-string] - file type of the data file, including regular expression notation (ex. '\*.dat')
* **myTitle** [array-string] - radargram title 
* **save_data** [logical] - controls whether to save (TRUE) or not save (FALSE) output .mat files
* **save_plots** [logical] - controls whether to save (TRUE) or not save (FALSE) output plots
* **filenames** [array-string] - input two strings corresponding to the file path location in which to save matched filter outputs and raw time domain data (inputted in that order)

##### Constants
* **startLat** [double] - starting latitude of transmitter
* **startLon** [double] - starting longitude of transmitter
* **SAMPLING_RATE** [integer] - sampling rate (in Hertz)
* **THICKNESS** [double] - ice thickness (in meters)
* **MAX_OFFSET** [double] - maximum antenna separation (in meters)
* **C** [double] - speed of light (in meters/second)

##### Loaded Files
The files in this section are already automatically loaded into the script, so no adjustments should be necessary. Only update the `load` function if the names change.
* **processing_inputs.mat** [struct] - contains the input parameters for the Greenland experiment used as an example
* **ref_chirp.mat** [array complex-double] - contains the reference chirp

#### Output
##### Files
* **allMatchFiltChirps.mat** [matrix complex-double] - a matrix containing the magnitudes of the matched filtered chirps at each antenna separation 
* **antennaSeparation.mat** [array-double] - an array containing the antenna separation distances (in meters) for each chirp in **allMatchFiltChirps**
* **gpsPts.mat** [array-double] - an array containing the GPS coordinate locations extracted from the last 100 bytes of each 8-s time-domain file  
##### Summary Plots
* *Radargram of all matched filtered chirps (x-axis = distance, y-axis = time)*
* *2D path plot of GPS position over time*
* *Radargram of all matched filtered chirps (x-axis = chirp index, y-axis = time)*

#### Internal Functions
`matchedFilterAlign` - a function that matched filters each received time-domain signal with the reference chirp, and applies coarse time-domain alignment to it.
* **Inputs**
  * **ref_chirp** [array complex-double] - the `ref_chirp.mat` file in `tools` containing the reference transmit signal
  * **fs** [integer] - the sampling rate (in Hertz)
  * **directory** [array-string] - contains the file path to the *mth* current time-domain data file in `data` 
  * **windowSize** [integer] - an internal variable that sets the plot viewing window
  * **dataType** - *see above*
  * **display** - *see above*
  * **gpsYN** - *see above*
  * **filenames** - *see above*
* **Outputs**
  * **matchFiltSub** [matrix complex-double] - a matrix containing the matched filtered, time-aligned chirps for all valid chirps in the current **directory** 
  * **numChirpsDetected** [integer] - the number of valid chirps detected in the current **directory**
  * **noChirpsFiles** [array-string] - set to **directory** if no valid chirps are detected
  * **chirpTimeFromEnd** [integer] - time difference between the direct path peak and the end of the recording (in seconds)
  * **gpsData** [struct] - contains the latitude, longitude, and time-stamp of the gps/transmitter extracted from the last 100 bytes of the SDR data
  * **noisePW** [array-double] - contains the noise power from the time-domain signal file (we subtract out the power of the direct path peak)
  * **noiseFlr** [array-double] - contains the average power (units of *V sqrt(s)*) from the time domain signal file
 
`readE312Data` - a function that parses the gps position data from the last 100 bytes of the time-domain data file
* **Inputs**
  * **directory** - *see above*
  * **dataType** - *see above*
  * **gpsYN** - *see above*
* **Outputs**
  * **data** [array-string] - is set equal to the value of **directory**, contains the SDR data filename 
  * **gpsData** - *see above*
. 
### Phase 2

<p align="center">
  <img width="675" height="350" src="https://github.com/bienert/Bistatic_Radar/blob/main/Misc/Images/phase2.png">
</p>

#### Description
The outputs from Phase 1 are loaded into Phase 2 as inputs. We first resample each chirp to a new direct path peak in the frequency domain (found using quadratic least squares peak estimation) in `parabolicTimeShift` (as shown in the first plot in the Phase 2 figure). Next, each chirp is upsampled and time-aligned to a template waveform. We form a group of chirps that reflected from the ice-sheet bed in the same approximate Fresnel zone as the current/center chirp. All the chirps in the Fresnel bin are then fine-scale phase-aligned to the center chirp (second plot in the Phase 2 figure). We use a normalized matched filter to compare each bin-chirp to the center chirp and remove chirps with poor correlation. Finally, we coherently sum and average all chirps in each Frensel bin, producing a single chirp with a visible bed echo to be used in further processing (third plot in the Phase 2 figure). 
#### Input
All input constants (in the **Constants** section) are automatically loaded from the `processing_inputs.mat` file to ensure the continuity of values throughout all phases of processing. See the [Variable Input Methods](#variable-input-methods) section for information on how to edit the file. All input controls (in the **Controls** section) should be entered manually.
##### Controls
* **filename** [array-string] - the filepath to an external folder to store generated plots
* **mod_freq** [integer] - the interval between indicies to save the central Fresnel bin chirp; this is useful for trouble-shooting
* **save_figs** [logical] - controls whether to save figures (TRUE) or not (FALSE). Useful for trouble-shooting; the final summary figures are saved regardless of the value of **save_figs**
* **addEnvelope** [logical] - determines whether to add apply an envelope function to the central Fresnel bin chirp (TRUE) in attempt to improve bin chirp cross-correlation, or not (FALSE). The default is **addEnvelope** = FALSE

##### Constants
* **SIGNAL_FREQUENCY** [integer] - the frequency of the transmit chirp (in Hertz)
* **ICE_THICKNESS** [double] - the thickness of the ice sheet (in meters)
* **ER_ICE** [double] - the real permitivity of the ice sheet
* **C** [double] - speed of light (in meters/second)
* **SAMPLING_RATE** [integer] - sampling rate (in Hertz)
* **UPSAMPLING_FACTOR** [integer] - the upsampling factor used in sinc interpolation upsampling
* **corr_scale_factor** [integer] - scales up the normalized matched filter outputs to more manageable values

##### Loaded Files
The files in this section are already automatically loaded into the script, so no adjustments should be necessary. Only update the `load` function if the names change.
* **allMatchFiltChirps** [matrix complex-double] - contains the matched-filtered signals outputted from *Phase 1*
* **antennaSeparation** [array-double] - the antenna separation distance (in meters) of each chirp in **allMatchFiltChirps**

#### Output
##### Files
* **coherentSum** [matrix complex-double] - matrix containing the coherently-summed chirp at each antenna separation
* **coherentAverage** [matrix complex-double] - matrix containing the coherently-averaged chirp (summed chirp / bin-length) at each antenna separation
* **numKeptChirps** [array-double] - an array containing the number of chirps kept from each Fresnel bin
##### Summary Plots
* *Number of kept chirps in each Fresnel bin vs. chirp index (also known as "antenna separation index")*
* *Number of discarded chirps in each Fresnel bin vs. chirp index*
* *Mean Fresnel bin correlation to center chirp vs. chirp index*
* *Mean Fresnel bin phase difference from center chirp vs. chirp index*
* *Radargram of chirp (x-axis = distance, y-axis = time)*

#### Internal Functions
`gauss` - a function used to create a guassian curve based on an inputed domain, standard deviation, and mean
* **Inputs**
  * **X** [array-double] - a time-domain range in seconds used to create a template waveform for center-chirp time alignment
  * **S** [double] - a standard deviation used to create the template waveform
  * **M** [double] - the mean used to create the template waveform
* **Outputs**
  * **Y** [array-double] - outputs a gaussian curve centered at **M**

`parabolicTimeShift` - a function that applies quadratic least squares peak estimation to the direct-path peak of the matched filtered output in **allMatchFiltChirps**, and resamples the waveform to this improved peak measurement
* **Inputs**
  * **SAMPLING_RATE** - *see above*
  * **allMatchFiltChirps** - *see above*
* **Outputs**
  * **allMatchFiltChirps_shifted** [matrix complex-double] - a matrix containing the time-shifted matched filter outputs from *Phase 1* at each antenna separation      

### Phase 3

<p align="center">
  <img width="675" height="350" src="https://github.com/bienert/Bistatic_Radar/blob/main/Misc/Images/phase3.png">
</p>

#### Description
Phase 3 first converts each chirp in **coherentAverage** to power (in Watts) as shown in the first diagram of the Phase 3 figure. It then calls `getPower` on each chirp, which uses quadratic least squares peak estimation to better approximate the bed echo peak. After resampling the waveform to this improved measurement (second diagram of the Phase 3 figure), the function calls `checkBasalPeak` to determine whether the basal echo peak has a time-delay and SNR that falls within user-defined thresholds. The ice-depth at the location of all valid chirps is then calculated. The final results are saved in **powerChirps** and **ice_depths**.  
#### Input
All input constants (in the **Constants** section) are automatically loaded from the `processing_inputs.mat` file to ensure the continuity of values throughout all phases of processing. See the [Variable Input Methods](#variable-input-methods) section for information on how to edit the file. All input controls (in the **Controls** section) should be entered manually.
##### Controls
* **error_check** [logical] - determines whether to include (TRUE) an error-checking algorithm that removes bed picks appearing at a depth discontinuous from surrounding data point values, or not (FALSE). See the description for the function `checkBasalPeak` for more information 
* **TIME_THRESHOLD** [double] - the maximum time difference between the previous basal pick and current one (in seconds)
* **NOISE_THRESHOLD** [double] - the minimum power intensity (in decibels) that a basal peak must be above to be considered valid and not noise
* **display** [logical] - controls whether a live basal-picking visualization is shown (TRUE) or not (FALSE). Does not affect the final output figures
##### Constants
* **ICE_THICKNESS** [double] - the thickness of the ice (in meters)
* **SDR_RESISTANCE** [double] - the resistance in the SDR instrument (in ohms)
* **ER_ICE** [double] - the real permitivity of the ice sheet
* **C** [double] - the speed of light (in meters/second)
* **SAMPLING_RATE** [integer] - the original sampling rate (in Hertz). Note that this should be the same sampling rate used in *Phase 2*
* **UPSAMPLING_RATE** [integer] - the upsampling factor used in *Phase 2*

##### Loaded Files
The files in this section are already automatically loaded into the script, so no adjustments should be necessary. Only update the `load` function if the names change.
* **coherentAverage** [matrix complex-double] - contains the coherently averaged chirp at each antenna separation distance
* **ref_chirp** [array complex-double] - contains the reference chirp
* **antennaSeparation** [array-double] - the antenna separation distance (in meters) of each chirp in **allMatchFiltChirps**

#### Output
##### Files
* **powerChirps** [array-double] - an array containing the basal pick powers at each antenna separation
* **ice_depths** [array-double] - an array containing the estimated ice depths at each antenna separation
##### Summary Plots
* *Basal echo power vs. along path distance* - a plot showing how the basal echo power magnitude changes with along path distance (1/2*antenna separation)
* *Ice depth vs. along path distance*

#### Internal Functions
`checkBasalPeak` - a function that checks whether a given chirp has a viable basal echo peak by comparing the time at which it occured to the time-stamp of the previous viable basal echo peak, and also comparing the basal echo peak's signal-to-noise-ratio to a threshold. Basal echo peaks with a time-delay from the previous peak greater than a user-defined threshold, or that have a SNR lower than the user-defined threshold, are removed. The very first basal echo peak is error-checked by comparing its time-domain location to the predicted delay. Since the ApRes ice-depth data used in **ICE_THICKNESS** was taken at the same coordinate-location as the first antenna separation distance, we expect this predicted delay time to be accurate. This error-checking function can be turned off by setting the **error_check** variable to be FALSE.
* **Inputs**
  * **error_check** [logical] - controls whether error-checking is applied (TRUE) or not (FALSE)
  * **ind** [integer] - the index of chirp in the **coherentSum** matrix. This is also the index of antenna separation distance in the **antennaSeparation** array for the chirp
  * **prev** [integer] - the index (in samples) of the previous basal echo peak in the last valid chirp
  * **pos** [integer] - the index (in samples) of the current basal echo peak
  * **TIME_THRESHOLD** - *see above*
  * **delayBed** [double] - the predicted time-delay between the start of the chirp and the bed echo peak
  * **SAMPLING_RATE** [integer] - the modified sampling rate that includes the upsampling factor (in Hertz)
  * **NOISE_THRESHOLD** - *see above*
  * **snr** [double] - the estimated signal-to-noise-ratio of the basal echo peak to surrounding noise (the SNR is calculated for only the second half of the chirp)
* **Outputs**
  * **error** [logical] - a value of TRUE indicates that the chirp has an unviable basal echo peak and will be eliminated from the final outputs; a value of FALSE indicates the chirp basal peak is viable
  * **prev** [integer] - the index (in samples) of the checked chirp's basal echo peak to be used as the comparator for the next chirp in **coherentSum** 

`getPower` - a function that applies quadratic least squares peak estimation to the basal peak of all viable chirps. It then finds, outputs, and plots the basal echo power peak amplitude. In addition, it can estimate the ice thickness (in meters) at each antenna separation offset. `getPower` calls the function `checkBasalPeak` inside of it. Power and ice depths from chirps that are not valid are not recorded (disable error-checking to capture these values).
* **Inputs** 
  * **error** [logical] - contains the result of **error** outputed by `checkBasalPeak`. A value of TRUE indicates the chirp does not have a viable basal peak, while a value of FALSE indicates the opposite
  * **chirp** [array-double] - the magnitude of the current chirp 
  * **ind** - *see above*
  * **prev** - *see above*
  * **TIME_THRESHOLD** - *see above*
  * **SAMPLING_RATE** - *see above*
  * **NOISE_THRESHOLD** - *see above*
  * **V** [double] - the speed of the EM-wave in ice (in meters/second)
  * **C** - *see above*
  * **ER_ICE** - *see above*
  * **ICE_THICKNESS** - *see above*
  * **distance** [double] - the antenna separation distance (in meters) of the current chirp
  * **disp** [logical] - a value of TRUE will turn on the basal power picking visualization (a plot of the chirp with the selected basal echo peak highlighted)
* **Outputs**
  * **power** [double] - the detected basal echo peak power for valid chirps (in Watts)
  * **peak_ind** [integer] - the index which the basal echo peak was detected (in samples)
  * **depth** [double] - the estimated ice depth (in meters) at the current (valid) chirp's antenna separation distance 
  * **prev_ind** [integer] - the index of the basal echo peak index to be for the next chirp's error checking
  * **snr** [double] - the estimated signal-to-noise-ratio for the current (valid) chirp's basal echo peak
  * **error** [logical] - returns whether the current chirp was viable (FALSE) or not (TRUE)

### Example
Check out our example/demo code in the **Example** folder. Be sure to check out the ReadMe that walks you through it.

### Cluster Computing
We advice processing large volumes of data on high-performance computing clusters. Please keep in mind, that we have not successfully been able to parallelize this processing chain. If you do run this, we recomend that you comment out most, if not all, of the internal plots generated in Phase 2. Phase 1 only has plots generated in **allMatchFiltChirps** which can be toggled on and off with *display*. 

### Questions
If you have questions or comments about our work, feel free to reach out to [dustin.M.Schroeder@stanford.edu](mailto:dustin.M.Schroeder@stanford.edu)!

 


