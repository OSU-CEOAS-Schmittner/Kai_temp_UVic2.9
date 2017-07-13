#!/bin/csh
#$ -e /home/kai/UVic2.9/MOBI1.9/nobio
#$ -o /home/kai/UVic2.9/MOBI1.9/nobio
cd /home/kai/UVic2.9/default_comb2/nobio
source /share/apps/lf6481/csh_laheyfort_setup
time ./UVic_ESCM > pr

