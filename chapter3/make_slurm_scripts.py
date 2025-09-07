import numpy as np
import glob
import os 
from datetime import datetime



samps = 1000
betas = []
files= []
alphas = [[.04]]
taus = [[125]]
time_str=datetime.now()
flogstr ="submit_log_1"+".log"
flog = open(flogstr,"a")
flog.write('date/time:\n')
flog.write(str(time_str)+'\n')
flog.write('taus:\n')
flog.write(str(taus)+'\n')
flog.write('alphas:\n')
flog.write(str(alphas)+'\n')
flog.write('Bacs.:\n')
flog.write(str('Bs')+'\n')
tau = 125
alpha = .04
betas = np.array([20,50,75,100,125,150,175,200])

for beta in betas:
	param_str="beta" +str(beta)
	f = open("run_het_sims_sde_"+param_str+".slurm","a")
	files.append("run_het_sims_sde_"+param_str+".slurm")
	f.write('#!/bin/bash\n')
	f.write('#SBATCH -J run_sims_Nb100'+param_str+'\n')
	f.write('#SBATCH -A FUSCO-SL3-CPU\n')
	f.write('#SBATCH --nodes=1\n')
	f.write('#SBATCH --ntasks=1\n')
	f.write('#SBATCH --cpus-per-task=1\n')


	f.write('#SBATCH --time=04:00:00\n')

	f.write('#SBATCH --mem=5980mb\n')
	f.write('#SBATCH --array=1-'+str(samps)+'\n')
	f.write('#SBATCH -p skylake\n')
	f.write('./etc/profile.d/modules.sh\n')
	f.write('module purge \n')
	f.write('module load rhel7/default-peta4\n')
	f.write('echo "This is job" $SLURM_ARRAY_TASK_ID\n')
	#f.write('g++ -o outanc_UDMp_'+param_str+'_$SLURM_ARRAY_TASK_ID phage_inf_coarse_het_multi.cpp -std=c++11\n')
	f.write('./outanc'+ ' -t '+str(tau)+' -a '+str(alpha)+ ' -b ' + str(beta)+' -i $SLURM_ARRAY_TASK_ID\n')
	f.close()

for f in files:
	beta = f.split('beta')[1]
	print(beta)
	flog.write("beta: "+ str(beta) +'\n')
	os.system('sbatch '+str(f)+ ' >> '+flogstr)
os.system('mkdir slurm_scripts\n')
os.system('mv *.slurm slurm_scripts/\n')

		


