# 2024/03/04: Michel Bechtold
# To compile lis, run the following steps
cd LISF/lis
# load the modules (right ones for each system will be automatically selected)
source KUL_LIS_modules
# run the configure and accept all defaults or change to 'debug' compilation if wanted
./configure 
# submit either the hortense or genius job scripts
# hortense (tier 1), need to swap to milan
module swap cluster/dodrio/cpu_milan
qsub compile_lis_hortense.pbs
# genius (tier 2)
sbatch compile_lis_genius.slurm
