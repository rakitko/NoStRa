NoStRa: 

1. In folder Generator:

a) edit config.txt
b) execute  generate_fields.R

The output will be:

If in config.txt,
compute_field_true_covs=1,
then the output file model.RData will be created in 
Generator/Out.

If in config.txt,
perform_enkf=1, 
then the  output file filter.RData will be created in 
Generator/Out.

If in config.txt,
compute_field_true_covs=-1  and  perform_enkf=-1, 
then the  output file fields_long.RData will be created in 
Generator/Out.
This option is intended to be used with M=1 to generate
a realization of the field xi for a long time (e.g. 60000 time steps).


2. Move the file Generator/Out/fields_long.RData to folder 
Statistics/Fields_long.
Then, in folder Statistics/Fields_long, 
execute the script Exam_fields_long.R.
-> Plots with time asnd space xi correlation functions will be created.
   kurtosis(xi) will be evaluated and written to the file out_xi_all.txt.


3. Move the file Generator/Out/model.RData to folder 
Statistics/NstatioFieldCovs.
Then, in folder Statistics/NstatioFieldCovs, 
execute the script ExamDSMoutput.R.
-> Plots with the secondary fields and xi will be created 
   as well as plots with randomly chosen spatial field correlations.

4. Move the file Generator/Out/filter.RData to folder 
Statistics/Statistics/EnKF-covs.
Then, in folder EnKF-covs, execute the script ExamEnKFcovs.R.
-> Various plots related to the EnKF's prior covariances will be created.


17 Aug 2017
A.Rakitko
M.Tsyrulnikov





