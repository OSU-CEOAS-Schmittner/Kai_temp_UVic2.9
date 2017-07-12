#!/bin/csh
# BE SURE To UPDATE YOUR WORKING DIRECTORY BELOW
#$ -e /home/changeme/UVic2.9/MOBI1.9
#$ -o /home/changeme/UVic2.9/MOBI1.9
cd /home/changeme/UVic2.9/MOBI1.9
source /share/apps/lf6481/csh_laheyfort_setup
time ./UVic_ESCM > pr

