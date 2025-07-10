1. Inside folder 'Scripts_Assemble_GPS_Data' run a script:
This load files from Data_GPS and saves assembled GPS.json into Data_2022-11-??__I24_Base for each day.

2. Then generate the full or slimmed integrated data in the respective folder. This puts GPS data, I-24 data, and energy models together to make an integrated data set.
"Scripts/generate_data_mvt_full" or "Scripts/generate_data_mvt_slim." This is saved into the folder Data/2022-11-??_MVT_Slim (or Full).

In that process, the script generate_mvt_data.m uses the energy model functions in the folder 'Models.' 

3. Then generate the macro fields (Fig. 3) and plot them in two steps from scripts in the folder 'Scripts_macroscopic_fields.'
	a. Run 'generate_macroscopic_fields_motion.m'
	b. Run 'plot_macroscopic_fields_motion.m'
	
4. Then one can generate Fig. 2, and Fig. SM2. using scripts in the folder 'Scripts/generate_data_samples.' This requires two steps.
	a. Run 'Scripts/generate_data_samples.m' first using the slim or full data. It produces mat-files (one per day) in a different Folder called 'Data_Analysis,' which contains 3 mat-files, one for each day. 
	b. Run 'Scripts/plot_AV_analysis.m' next to generate the results from the article (figure 2, figure SM2, and Figure SM3). 










	
 




