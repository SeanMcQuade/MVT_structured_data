
These scripts and data are from the CIRCLES consortium MegaVanderTest experiment that was conducted in Nashville on highway I-24 during the week of November 14th 2022. 
This repository contains data recorded from the partially automated vehicles, called the GPS data (this includes GPS location data recoreded by the 
installed raspberri pi and several flags to indicate the state of the vehicle) and data recorded by the I-24 MOTION observatory.

# Contents
- [Generate integrated data.](#tag1)
- [Plot data and results.](#tag2)
- [Websites](#tag3)

## Scripts to generate the integrated data set.

Step 0: Have the correct folder structure before you begin. The main folder '/' must create the following subfolders /Data, /Models, and /Scripts. 
There should be 10 folders in the Data/ folder.
Three correspond to 11/16/2022: "Data_2022-11-16__I24_Base", "Data_2022-11-16__I24_Slim", "Data_2022-11-16__I24_Full",
Three correspond to 11/17/2022: "Data_2022-11-17__I24_Base", "Data_2022-11-17__I24_Slim", "Data_2022-11-17__I24_Full", and
Three correspond to 11/18/2022: "Data_2022-11-18__I24_Base", "Data_2022-11-18__I24_Slim", "Data_2022-11-18__I24_Full".
The tenth folder is called "Data_GPS".

The base data .json files go into "Data_2022-11-1*__I24_Base" then create the folders "Data_2022-11-16__I24_Slim" and "Data_2022-11-16__I24_Full". The slim and full data writes to these folders respectively.

Step 1: Run "Scripts\generate_data_mvt_full" or "Scripts\generate_data_mvt_slim.":This load files from Data_GPS and saves assembled GPS.json into Data_2022-11-??__I24_Base for each day. 

Step 2: Once the slimmed or full data is generated, run "Scripts\generate_data_samples.m." It produces mat-files in the Folder called 'Data_Analysis,' one for each day.

Step 3: Run 'Scripts\generate_macroscopic_fields.m.' This will save macoscopic fields data to the "Data\Data_for_Figures" folder as a matlab file ".mat."


## Scripts to plot figures from the article.
Step 4: Run "Scripts\plot_AV_stats.m' to generate the results from the article (figure 2, figure SM2, and Figure SM3).

Step 5: Run "Scripts\plot_macroscopic_stats.m' to generate the macroscopic fields figures from the article (figure 3, and SM 5, as well as additional fields).

## Websites
[Visit the CIRCLES consortium website](https://circles-consortium.github.io/)

[This repository is licensed under the BSD3 license](https://opensource.org/license/bsd-3-clause)
