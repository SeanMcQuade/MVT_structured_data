
These scripts and data are from the CIRCLES consortium MegaVanderTest experiment that was conducted in Nashville on highway I-24 during the week of November 14th 2022. 
This repository contains data recorded from the partially automated vehicles, called the GPS data (this includes GPS location data recoreded by the 
installed raspberri pi and several flags to indicate the state of the vehicle) and data recorded by the I-24 MOTION observatory.

# Contents
- [Generate integrated data.](#tag1)
- [Plot data and results.](#tag2)
- [Websites](#tag3)

## Scripts to generate the integrated data set.

### Step 0: Correct folder structure before you begin. 

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


### Step 1: run some scripts to generate the data

You have a couple options here.

---
Run "Scripts\generate_data_mvt_full" or "Scripts\generate_data_mvt_slim.":This load files from `../data/cars/gps` and saves assembled GPS.json into Data_2022-11-??__I24_Base for each day. 

Step 2: Once the slimmed or full data is generated, run "Scripts\generate_data_samples.m." It produces mat-files in the Folder called 'Data_Analysis,' one for each day.

Step 3: Run 'Scripts\generate_macroscopic_fields.m.' This will save macoscopic fields data to the "Data\Data_for_Figures" folder as a matlab file ".mat."


## Scripts to plot figures from the article.
Step 4: Run "Scripts\plot_AV_stats.m' to generate the results from the article (figure 2, figure SM2, and Figure SM3).

Step 5: Run "Scripts\plot_macroscopic_stats.m' to generate the macroscopic fields figures from the article (figure 3, and SM 5, as well as additional fields).

## Websites
[Visit the CIRCLES consortium website](https://circles-consortium.github.io/)

[This repository is licensed under the BSD3 license](https://opensource.org/license/bsd-3-clause)
