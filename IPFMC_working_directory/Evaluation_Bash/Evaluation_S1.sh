taskset -c 1 nohup python ../Data_Preprocess/IPFMC_strategy_1/Prepare_SMs_SC3ARIS_CNV.py > CNout.log &
taskset -c 2 nohup python ../Data_Preprocess/IPFMC_strategy_1/Prepare_SMs_SC3ARIS_mRNA.py > mRout.log &
taskset -c 3 nohup python ../Data_Preprocess/IPFMC_strategy_1/Prepare_SMs_SC3ARIS_miRNA.py > miRout.log &
taskset -c 4 nohup python ../Data_Preprocess/IPFMC_strategy_1/Prepare_SMs_SC3ARIS_Methy.py > Metout.log &
