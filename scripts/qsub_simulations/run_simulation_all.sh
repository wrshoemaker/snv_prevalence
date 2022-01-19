#!/bin/bash




#declare -a all_f_mean=(0.001 0.00199526 0.00398107 0.00794328 0.01584893 0.03162278 0.06309573 0.12589254 0.25118864 0.50118723)
#declare -a all_beta=(3.16227766e-03 7.74263683e-03 1.89573565e-02 4.64158883e-02 1.13646367e-01 2.78255940e-01 6.81292069e-01 1.66810054e+00 4.08423865e+00 1.00000000e+01)

#for f_mean in "${all_f_mean[@]}"
#do
#  for beta in "${all_beta[@]}"
#  do
#    qsub /u/home/w/wrshoema/project-ngarud/snv_prevalence/scripts/qsub_simulations/qsub_simulation.sh ${f_mean} ${beta}
#  done
#done




#0.001 0.00138709 0.00192402 0.0026688  0.00370187 0.00513483 0.00712249 0.00987954 0.01370383 0.01900848 0.02636651 0.03657277 0.0507298  0.0703669  0.09760536 0.13538762 0.18779508 0.26048905 0.36132228 0.50118723


declare -a all_f_mean=(0.09760536 0.13538762 0.18779508 0.26048905 0.36132228 0.50118723)
declare -a all_beta=(1.00000000e-03 1.62377674e-03 2.63665090e-03 4.28133240e-03 6.95192796e-03 1.12883789e-02 1.83298071e-02 2.97635144e-02 4.83293024e-02 7.84759970e-02 1.27427499e-01 2.06913808e-01 3.35981829e-01 5.45559478e-01 8.85866790e-01 1.43844989e+00 2.33572147e+00 3.79269019e+00 6.15848211e+00 1.00000000e+01)



for f_mean in "${all_f_mean[@]}"
do
  for beta in "${all_beta[@]}"
  do
    qsub /u/home/w/wrshoema/project-ngarud/snv_prevalence/scripts/qsub_simulations/qsub_simulation.sh ${f_mean} ${beta}
  done
done




#/u/home/w/wrshoema/project-ngarud/snv_prevalence/scripts/qsub_prevalence/qsub_predicted_prevalence.sh 0.1 1.0
