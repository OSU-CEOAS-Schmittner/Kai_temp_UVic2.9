#!/bin/csh
#$ -e /raid24/aschmitt/UVic2.9/MOBI1.9
#$ -o /raid24/aschmitt/UVic2.9/MOBI1.9
cd /raid24/aschmitt/UVic2.9/MOBI1.9
source /share/apps/lf6481/csh_laheyfort_setup
time ./UVic_ESCM > pr

