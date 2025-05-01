In this folder you will find folders to generate several data sets. They should be run in the order presented:
1. ./Generate_integrated_data_slim/generate_data_mvt_slim.m
	This script generates the data from the slim data from GPS .json file and the I-24 MOTION .json files. They are written to the folder ./2022-11-1*__MVT_Data_Slim
2. ./Generate_integrated_data_slim/generate_data_samples.m from data in ./2022-11-1*__MVT_Data_Slim They are written to the same folder (./2022-11-1*__MVT_Data_Slim)
	This script generates integraged data samples used to produce the main results figure of the article (Fig. 2).
3. ./Generate_macroscopic_fields/generate_macroscopic_fields_motion.m
	This script generates macroscopic fields from the output files from data in ./2022-11-1*__MVT_Data_Slim.

Also in this folder is plotting scripts to be used after the data has been generated:
4. /Plot_FC_results/Plot_AV_stats_effective_FC.m uses data samples in ./2022-11-1*__MVT_Data_Slim to produce the main results figure of the article (Fig. 2).
5. /Plot_FC_results/Plot_AV_stats.m uses data samples in ./2022-11-1*__MVT_Data_Slim to produce the supplemental results figure of the article (Fig. SM 2).
6. /Plot_macroscopic_fields_slim/mvt_v21_plot_macroscopic_fields_motion_nature requires data in ./2022-11-1*__MVT_Data_Slim to produce the time-space diagram with AV trajectories overlaid (figure 3) from the article.








	
 




