taskset -c 1 nohup python ../Downsteam_Analysis/Spec_Surv_Complete.py > Survival.log &
taskset -c 2 nohup python ../Downsteam_Analysis/TrueLabel_Complete.py > Truelabel.log &