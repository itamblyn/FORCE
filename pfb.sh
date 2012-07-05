#!/bin/tcsh

tail -n +11 TRAJEC.cbn > tmp.cbn

wc -l tmp.cbn FTRAJECTORY

paste FTRAJECTORY tmp.cbn > tmp.tmp

awk '{print $2, $3, $4, $8, $9, $10, $15}' tmp.tmp > position_force_bonding

rm tmp.cbn tmp.tmp
