#!/bin/bash
#####################################################
if [ -z "$1" ] ;then
    echo "No argument supplied"
    echo enter with the phase: 1 to 5
    exit
fi
#####################################################
fase=${1}

mesh1=x1.10242    # 240km 
mesh2=x1.40962    # 120km
mesh3=x1.163842   #  60km
mesh4=x1.655362   #  30km 
mesh5=x1.2621442  #  15km

#new names
exp6=GFdef 
#####################################################
mesh_target=${mesh1}


CPORG=0
GUSTF=0 # sea spray make this 2
SUB3D=0

expname=shocmf

run=6
nruns=6
while [ $run -le $nruns ]
do
#----
if [ $run == 0  ]; then exp=${exp0} ; mesh=${mesh_target}; fi
if [ $run == 1  ]; then exp=${exp1} ; mesh=${mesh_target}; fi
if [ $run == 2  ]; then exp=${exp2} ; mesh=${mesh_target}; fi
if [ $run == 3  ]; then exp=${exp3} ; mesh=${mesh_target}; fi
if [ $run == 4  ]; then exp=${exp4} ; mesh=${mesh_target}; fi
if [ $run == 5  ]; then exp=${exp5} ; mesh=${mesh_target}; fi
if [ $run == 6  ]; then exp=${exp6} ; mesh=${mesh_target}; fi
if [ $run == 7  ]; then exp=${exp7} ; mesh=${mesh_target}; fi
if [ $run == 8  ]; then exp=${exp8} ; mesh=${mesh_target}; fi
if [ $run == 9  ]; then exp=${exp9} ; mesh=${mesh_target}; fi
if [ $run == 10 ]; then exp=${exp10}; mesh=${mesh_target}; fi
if [ $run == 11 ]; then exp=${exp11}; mesh=${mesh_target}; fi

make_gf_nml.bash ${exp} 
run.bash ${fase} ${exp} ${mesh} ${expname} ${CPORG} ${GUSTF} ${SUB3D} > ${exp}_${mesh}_${expname}.out & 
echo 'run with exp/fase/mesh=' ${exp} ${fase} ${mesh} ${expname}
sleep 5

#----
run=`expr $run + 2`
done

/bin/cp exp.bash  ${mesh_target}_${expname}_${exp}_exp.bash
chmod +x ${mesh_target}_${expname}_${exp}_exp.bash

exit





