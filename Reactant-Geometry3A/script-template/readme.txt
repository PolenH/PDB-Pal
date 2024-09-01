For this script all you need to do is call script.sh. Command will look like bash script.sh -nc -2 for instance, then call sbatch mm_orca_job.sh.
The -nc flag is for the given charge of the ligand, only thing that has to be manually inputed. This script
still has issues with ligands that have structures that need to remain planar due to Orca. Next will be checking the file manually
this will let us know if we need to relax more atoms which requires more manual input. This needs to be workshopped.