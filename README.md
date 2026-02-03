
These scripts and data are from the CIRCLES consortium MegaVanderTest experiment that was conducted in Nashville on highway I-24 during the week of November 14th 2022.  This repository contains data recorded from the partially automated vehicles, called the GPS data (this includes GPS location data recoreded by the installed raspberri pi and several flags to indicate the state of the vehicle) and data recorded by the I-24 MOTION observatory.

These scripts provide analysis regarding the energy usage of vehicles during the CIRCLES consortium MegaVanderTest field experiments. The field experiments were conducted in Nashville on highway I-24 during the week of November 14th 2022. 

The repository refers to data recorded from the I-24 MOTION observatory, as well as the CIRCLES partially automated vehicles. Those data must be obtained elsewhere, and will be privately stored or available only to CIRCLES Team members until the paper is released.

For the peer review purposes, use branch: peer-review (https://github.com/SeanMcQuade/MVT_structured_data/tree/peer-review)

# Contents
- [Learn about the CIRCLES data.](#about-circles-data-analyzed-with-this-software)
- [Data installation](#data-install)
- [Plot data and results.](#plot-results)
- [(Advanced) Bootrstrap: Generate integrated data.](#advanced-bootstrap)
- [Websites](#tag3)

### 1. System Requirements

To carry out the analysis, the following software requirements are needed:

#### Software Dependencies:

- Install [MATLAB](https://mathworks.com/)
- Only `Matlab` must be installed as a required toolbox
- MATLAB v2025a and v2025b have both been tested
- Tested on Mac Sequoia 15.6, Windows 11

#### Installation Guide
- Install MATLAB (approximately 20 minutes, depending on download speeds)
- Fetch the data (will be made publicly available upon publication). The extraction should be to a hard drive with at least 700GB of available space, to ensure enough space to generate additional files.

#### Demo

##### Instructions to run the demo:

- Navigate to the folder `MVT_structured_data/Scripts`.
- Run the file `run_all_scripts.m`

The execution of this file will take significant time. It will reproduce the style of plot for the data provided. 

- Generate structured collections of data samples for this day from the `results/` folder
- Generate the macroscopic fields (from the provided data)
- Plot the macroscopic fields (from the provided data) with the GPS data from the cars (across the entire day)
- Plot the microscopic fields of all cars (from above)
- Carry out the analysis of fuel usage and comparisons (from the provided data)
- Plot the fuel usage and comparisons (from above)

A subset of these plots are synthesized in `results/figures/2022-11-17` in the released data set (available upon conclusion of peer review).

##### Expected output: 

Resulting files should be generated in the following folders for day 2022-11-17

###### Found in `results/figures/2022-11-17`
```
fields_motion_2022-11-17.mat
fig_2_fuel_results_effective_645_915.fig
fig_2_fuel_results_effective_645_915.png
fig_3_fuel_results_effective_mean_median_645_915.fig
fig_3_fuel_results_effective_mean_median_645_915.png
fig_field_20221117_west_laneall_motion_F_av_nature_large.png
fig_field_20221117_west_laneall_motion_Phi_av_nature_large.png
fig_field_20221117_west_laneall_motion_Psi_av_nature_large.png
fig_field_20221117_west_laneall_motion_Q_av_nature_large.png
fig_field_20221117_west_laneall_motion_Rho_av_nature_large.png
fig_field_20221117_west_laneall_motion_U_av_nature_large.png
fig_motion_trajectories_20221117_west_laneall_lowres.png
fig_motion_trajectories_20221117_west_laneall_zoom_lowres.png
fig_motion_trajectories_20221117_west_laneall_zoomwin_lowres.png
fig_SM2_vehicle_samples_counts_effective_mean_median_645_915.fig
fig_SM2_vehicle_samples_counts_effective_mean_median_645_915.png
samples_for_distance_analysis_17.mat
```

An example image that should be produced represents the macroscopic flow rates, with overlay  is `results/figures/2022-11-17/fig_field_20221117_west_laneall_motion_Psi_av_nature_large.png` which describes the bulk fuel consumption with overlay of GPS data from our control cars, indicating when their control was active (or not) during their drives.

![Bulk Fuel results (partial) with only a subset of data from 2022-11-17](../results/figures/2022-11-17/fig_field_20221117_west_laneall_motion_Psi_av_nature_large.png)

**Note** there may be minor errors or warnings thrown, since the data pipeline is intended to reproduce exact figures with comparable max/min values and colors across multiple plots. If only one day with a subset of data into that anaysis, is included, plots may have

- what seems to be missing large portions on the left/right of presented data
- what seems to indicate that axes are 'zoomed out'

These are artifacts of axes bounds that are normalized across multiple days for comparison, and thus the appropriate approach is to open those .fig files as plots with MATLAB, and then zoom in to explore.

###### Found in `results/slim/2022-11-17`

This will include a single .mat file that represents binary data in MATLAB format for quick review and analysis in subsequent plots.

```
I-24MOTION_2022-11-17_07-59-59_reduced.mat
```

###### Notes on warning messages

*Note* several output messages will show in the MATLAB window that show warnings for additional legend entries that are not used. This is due to only a subset of data being shown for the demonstration data.

```
Warning: Ignoring extra legend entries. 
> In legend>process_inputs (line 575)
In legend>make_legend (line 294)
In legend (line 245)
In plot_AV_analysis>plot_one_sided (line 473)
In plot_AV_analysis (line 310)
In reproduce_plots (line 27) 
```

##### Expected Runtime

It should take approximately 4-5 minutes or faster to regenerate plots if all intermediate data already are downloaded, and the `slim/` folder exists.


## About CIRCLES data analyzed with this software

### About data from I-24 MOTION

The `results/slim` data include trajectories recorded and processed by the I-24 MOTION observatory. Primarily the data from I-24 MOTION are comprised of position and speed trajectories from all vehicles detected on the roadway by that observatory. The `results/slim` data align the location and state of each of the CIRCLES control cars includes GPS location data as well as data from on-board vehicle informatics recorded by custom hardware installed with a Raspberri Pi. Those data include information regarding to indicate the control and system state of the vehicle. 

### About data from GPS and the CIRCLES Cars
Two kinds of data were collected from the CIRCLES Cars, and found in the `results/gps` folder in the data archive. Each day provides a separate standalone file for the cars that ran on that day. Information includes GPS data, aligned with timeseries information from data collected directly from the car.

#### Data collection from CIRCLES CARS with team-installed GPS Sensors
Information from team-installed GPS sensors was collected at 10-Hz. These data include the position and speed of the vehicle at each sample point.

#### Data Collection from CIRCLES Cars with team-designed on-board data collection
Information from team-installed computers that interface with the Controller Area Network (CAN) were critical to sensing and control of the experiment cars. These data are aligned with the raw GPS information to provide the state of the vehicle at that time (speed, assigned lane of travel, desired cruise control set point, etc.). 

## Data Install

Download the data to your computer, it will be called either `data` or `results`. In the same folder that contains the data or results, clone this repository.

## Reproduce Plots: Scripts to plot figures from the article.

Ensure you have downloaded all the data, and stored it according to the structure described in the next step. 

### Step 0: Correct folder structure for reproducing plots. 

This git repository should be a *sibling* folder to the `results` folder. 

```
cd MVT_structured_data
ls ..
```
This command should show you 

```
results/
MVT_structured_data
```

You may also have folders like `data` if you have downloaded the bootstrapping/base data.

If you want to confirm what it should look like from the 'parent' directory:

```
[Parent]
  - MVT_structured_data
  - results/
  | - slim
    | - 2022-11-16
      | - I-24MOTION_2022-11-16_05-59-59.json
      | - I-24MOTION_2022-11-16_06-09-59.json
      | - ...
    | - 2022-11-17
      | - ...
    | - 2022-11-18
      | - ...
  | - gps
    | - CIRCLES_GPS_10Hz_2022-11-16.json
    | - CIRCLES_GPS_10Hz_2022-11-17.json
    | - CIRCLES_GPS_10Hz_2022-11-18.json
  | - README.md
```

Once you run the scripts, additional folders and files will be produced in the `results` folder.

### Step 1: Reproduce the plots
To reproduce all the plots for all days, simply run 

```
cd MVT_structured_data/Scripts
```
and in that folder in matlab, run

```
reproduce_plots.m
```

This will run the below scripts for all days as well.

Running `Scripts\plot_AV_analysis.m` generates the results from the article (figure 2, figure SM2, and Figure SM3).

Running``"Scripts\plot_macroscopic_fields.m` to generate the macroscopic fields figures from the article (figure 3, and SM 5, as well as additional fields).

### Step 2: Examine figures

The outputs in `results/figures` provide reproductions of the figures used in the main graphics in the paper.

## Advanced Bootstrap
How to generate the integrated data set. (advanced only)

***Note*** This step requires a different dataset to begin, and generates the `slim` and `full` results that are part of the release. These steps perform alignment of vehicle and I-24 MOTION data, from base files from I-24 MOTION and original GPS and vehicle CAN data files that are assigned by each car.

This step is ADVANCED and OPTIONAL and is included mainly to provide the algorithmic insights to anyone interested in how those data are produced and aligned.

### Step 0 (Advanced): Correct folder structure for bootstrap data synthesis. 

This git repository should be a *sibling* folder to the data folder. 

```
cd MVT_structured_data
ls ..
```
This command should show you 

```
data/
MVT_structured_data
```

A folder named `results` will be created as part of these scripts.

If you want to confirm what it should look like from the 'parent' directory:

```
[Parent]
  - data/
  | - cars
    | - cars_gps
      | - circles_v2_1_car1.csv
      | - circles_v2_1_car2.csv
      | - circles_v2_1_car3.csv
      ...
    | - cars_vins.csv
    | - veh_ping_ ...
    | - ...
  | - i24motion
    | - 2022-11-16
      | - 64888dc ....wed_0_00.json
      | - 64888dc ....wed_0_01.json
      | - 64888dc ....wed_0_02.json
      | - ...
    | - 2022-11-17
      | - ...
    | - 2022-11-18
      | - ...
  | - README.md
  - MVT_structured_data
  - results
```

### Step 1 (Advanced): run scripts to generate the data

Run `Scripts/run_all_scripts.m` which begins with the `base` data and produces all intermediate files, aligns GPS information, and synthesizes intermediate storage in `.mat` format that can ease additional plotting and analysis by future researchers.

Some information on each of these scripts is below.

#### `Scripts/assemble_data_GPS.m`

Generates files in `results/gps` named `CIRCLES_GPS_10Hz_2022-11-{16,17,18}.json`, which combine individual vehicle data from each car (found in `data/cars` and `data/cars/cars_gps`).

#### `Scripts/generate_data_mvt_full.m`

This load files from `../data/cars/cars_gps` along with data from I-24 MOTION in `../data/i24motion/2022-11-{16,17,18}/*.json` and saves assembled GPS.json into `results/full/2022-11-{16,17,18}/I-24MOTION_2022-11-*.json` for each day. 

The "full" dataset includes trajectories from both westbound and eastbound lanes of I-24 Motion, though only westbound lanes are included in our analysis. Thus, there is no reason to utilize the `full` dataset for this paper, though the results may be of use to other researchers.


#### `Scripts\generate_data_mvt_slim.`

This load files from `../data/cars/cars_gps` along with data from I-24 MOTION in `../data/i24motion/2022-11-{16,17,18}/*.json` and saves assembled GPS.json into `results/slim/2022-11-{16,17,18}/I-24MOTION_2022-11-*.json` for each day. 

The "slim" dataset includes trajectories from *only* westbound lanes of I-24 Motion, and also removes several fields from the resulting datafile, which are not utilized in our analysis. The `slim` dataset is the preferred dataset for release, given its slightly smaller size.

#### `Scripts\generate_data_samples.m`

Once the slimmed or full data is generated, this script produces MATLAB `.mat` files in the `full` (or `slim`) folder for each day. These data samples include the speed, position, vehicle class, fuel rates, etc. that are needed to compute the macroscopic fuel usage results.

#### `Scripts\generate_macroscopic_fields.m`

This scripts utilizes the `generate_data_samples` function to calculate macoscopic field data, preparing for plotting and AV influence analysis. The results are saved in `results/{slim,full}/2022-11-{16,17,18}/*.mat`. 

#### `Scripts/plot_macroscopic_fields.m`

This script produces Figure 3, and SM 5, as well as additional fields.


#### (Optional) `Scripts/plot_microscopic_trajectories.m`

***Note*** Plotting these trajectories may cause your computer to run slowly, or in some cases to crash, if there are challenges with memory or data access. These are supplemental figures that show the microscopic trajectories for each day, and thus represent significant numbers of points and can render graphing functions inadequate.

#### `Scripts/plot_AV_analysis.m`

This script produces Figure 2, Figure SM2, and Figure SM3.

## Websites
[Visit the CIRCLES consortium website](https://circles-consortium.github.io/)

[This repository is licensed under the BSD3 license](https://opensource.org/license/bsd-3-clause)
