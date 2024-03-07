# 2024/03/04: Michel Bechtold
# To compile ldt, run the following steps
cd LISF/ldt
# load the modules (right ones for each system will be automatically selected)
source KUL_LDT_modules
# run the configure and accept all defaults or change to 'debug' compilation if wanted
./configure 
# submit either the hortense or genius job scripts
# hortense (tier 1), need to swap to milan
module swap cluster/dodrio/cpu_milan
qsub compile_ldt_hortense.pbs
# genius (tier 2)
sbatch compile_ldt_genius.slurm
