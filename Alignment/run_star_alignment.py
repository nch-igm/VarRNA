import sys, os



samples = ['sample_test']
for sample in samples:
    cmd = f'echo "bash star_alignment.sh {sample}" | qsub -pe smp 20 -cwd -j y'
    os.system(cmd)
