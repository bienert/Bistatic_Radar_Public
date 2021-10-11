# Processing Scripts Folder ReadMe
Please note that this repository accompanies a paper in review titled "Large Offset Bistatic Radar Sounding of Glaciers using Processing-Based Synchronization". Please cite this paper to comply with our license when redistributing transformed versions of this work (e.g. publishing a paper on experiments which which used some code from this repo). For example: N. Bienert et al., "Post-Processing Synchronized Bistatic Radar for Long Offset Glacial Sounding," IEEE TRANS. Geosci. Remote Sens., in review.

## Contents

### Full_Processing 
**GO HERE** This folder contains the scripts and instructions for transforming raw unsynchronized bistatic radar data into bed picks and attenuation estimates.

### Simple_Processing
This folder contains scripts for looking at time domain, frequency domain, and matched filtered data from the SDRs for troubleshooting and testing purposes when collecting data in the field. These scripts only accept a single file at a time and do not implement synchronization.

### functions
This folder contains functions used by scripts in both the *Full_Processing* and *Simple_Processing* folders.

### data 
Using our data: Download our [data](https://www.usap-dc.org/view/dataset/601472) into the *raw_data* folder, maintaining the file structure:
```
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

Your data: Place data you collected in the *raw_data* folder. You will need to adjust the *setup1_inputs_struct.m* in the *Full_Processing* folder to accomodate the data.
