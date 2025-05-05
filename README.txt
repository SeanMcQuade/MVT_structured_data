1. Inside folder 'Scripts_Assemble_GPS_Data' run a script:
This load files from Data_GPS and saves assembled GPS.json into Data_2022-11-??__I24_Base for each day.

2. Then generate the full or slimmed integrated data in the respective folder. This puts GPS data, I-24 data, and energy models together to make an integrated data set.
"Scripts_Integrate_Data_Full" or "Scripts_Integrate_Data_Slim"
in each, there is one script that does the following: (integrate_mvt_data.m) load files from Data_2022-11-??__I24_Base and saves into Data_2022-11-??_MVT_Slim (Full)
and it copies the GPS data file from Data_2022-11-??__I24_Base into Data_2022-11-??_MVT_Slim (Full)

In that process, the file integrate_mvt_data.m uses the energy model functions in the folder 'Models_energy.' The GPS data from step 1. is needed for plotting as well as for MVT data assembly in step 2.
Scripts_Macroscopic_Fields .m file inside loads data only from one directory 'Data_2022-11-??_MVT_Slim'

3. Then one can generate Fig. 2, and Fig. SM2. using scripts in the folder 'Scripts_Data_Analysis.' This requires two steps.
	a. Run 'generate_data_samples.m' first using the slim data. It produces mat-files (one per day) in a different Folder called 'Data_Analysis,' which contains 3 mat-files, one for each day. 
	b. Run 'Plot_AV_stats_effective_FC.m' next to generate the results from the article (figure 2). Run 'Plot_AV_FC_stats.m' next to generate the results for three different measures of Fuel Consumption from the article (figure SM2)n
	
4. Then one can generate the macro fields (Fig. 3) and plot them in two steps from scripts in the folder 'Scripts_macroscopic_fields.'
	a. Run 'generate_macroscopic_fields_motion.m'
	b. Run 'plot_macroscopic_fields_motion.m'











	
 




