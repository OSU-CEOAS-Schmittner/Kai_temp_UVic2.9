#!/bin/csh
#$ -e /home/kai/UVic2.9/default_comb2/nobio_cons_kgm
#$ -o /home/kai/UVic2.9/default_comb2/nobio_cons_kgm
cd /home/kai/UVic2.9/default_comb2/nobio_cons_kgm
source /share/apps/lf6481/csh_laheyfort_setup
time ./UVic_ESCM > pr

