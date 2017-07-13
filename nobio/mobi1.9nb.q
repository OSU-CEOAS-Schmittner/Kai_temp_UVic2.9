#!/bin/csh
#$ -e /raid24/aho/UVic2.9/MOBI1.9/nobio
#$ -o /raid24/aho/UVic2.9/MOBI1.9/nobio
cd /raid24/aho/UVic2.9/default_comb2/nobio
source /share/apps/lf6481/csh_laheyfort_setup
time ./UVic_ESCM > pr

