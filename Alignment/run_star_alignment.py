import sys, os


samples = ['IGMCH0018', 'IGMCH0131', 'IGMCH0140', 'IGMCH0143', 'IGMCH0159', 'IGMCH0201']

samples = ['IGMCH0018']
for sample in samples:
    cmd = f'echo "bash star_alignment.sh {sample}" | qsub -pe smp 20 -cwd -j y'
    os.system(cmd)
