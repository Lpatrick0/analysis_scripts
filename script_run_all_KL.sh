# have 2 folders, one with reference and another called target. Edit path to both. Set up folder for output files.

python KL_run_all_test.py orxA_psi/psi_hist_0.dat orzA_psi/psi_hist_0.dat output/psi_0.dat

rep=1
while [ $rep -lt 286 ];
do
        python KL_run_all_test.py orxA_psi/psi_hist_$rep.dat orzA_psi/psi_hist_$rep.dat output/psi_$rep.dat
	echo $rep
	rep=$((rep+1))

done
