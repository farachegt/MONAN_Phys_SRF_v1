#!/bin/bash
DIR0=/mnt/beegfs/guilherme.farache/runs/MONAN/model
DIR1=/mnt/beegfs/guilherme.farache/runs/MONAN/posprocess
DIR_GRID_SPEC=/mnt/beegfs/guilherme.farache/runs/MONAN/model/datain
DIR_CONVERT_MPAS=/mnt/beegfs/guilherme.farache/runs/convert_mpas_pkubota/convert_mpas_v0.1.0_egeon.gnu940/
local_shells=/mnt/beegfs/guilherme.farache/runs/scripts_shells
######################################################
#PARTITION AVAIL  TIMELIMIT  NODES  STATE NODELIST
#batch        up    4:00:00     33   idle n[01-33]
#proc         up    4:00:00      3   idle proc[01-03]
#OPER         up    3:00:00     15   idle n[01-15]
#PESQ1        up   12:00:00      8   idle n[16-23]
#PESQ2        up   12:00:00      8   idle n[24-31]
#PESQ3        up    2:00:00      2   idle n[32-33]
#DEBUG        up      15:00      1   idle proc01
#####################################################
if [ -z "$1" ] ;then
    echo "No argument supplied"
    echo enter with the phase: 1 to 5
    exit
fi
if [ -z "$2" ] ;then
    echo enter with GF configuration
    exit
fi
if [ -z "$3" ] ;then
    echo enter with mesh: x1.10242, x1.40962, x1.655362 _or_ x1.2621442
    exit
fi
if [ -z "$4" ] ;then
    echo enter with the expname
    exit
fi
if [ -z "$5" ] ;then
    echo enter with CPORG
    exit
fi
if [ -z "$6" ] ;then
    echo enter with GUSTF
    exit
fi

#####################################################
 
 JobName=SHOCMF
 Partition=batch
 RequestedTime=08:00:00
 output_interval_dc=1:00:00
 #Partition=PESQ3
 monthly_run=0
 change_gf_namelist=0
 #-- controls convertion of native MONAN/MPAS netcdf files
 doconvert=1

 fase=${1}
 exp=${4}_${2}
 

 if [ $fase == 1 ] ; then 
   exec=init_atmosphere_model
   #exec=stable_init_atmosphere_model
 fi
 if [ $fase == 2 ] ; then 
   exec=atmosphere_model
   #exec=stable_atmosphere_model
   #--- create the GF namelist
   make_gf_nml.bash ${2}
 fi
 
 generic_namelists=MONAN_generic_namelists
#generic_namelists=Stable_generic_namelists

mesh=${3}
#mesh=x1.10242    # 240km <<<<<<<<<<<
#mesh=x1.40962    # 120km
#mesh=x1.163842   #  60km
#mesh=x1.256002   #  48km 
#mesh=x1.655362   #  30km <<<<<<<<<<<
#mesh=x1.1024002  #  24km
#mesh=x1.2621442  #  15km

#para alta resolucao:
#config_physics_suite=convection_permitting_monan

#config_physics_suite=convection_permitting

#para baixa resolucao:
config_physics_suite=mesoscale_reference_monan
#config_physics_suite=mesoscale_reference
config_dt_cu=00:15:00
config_radt_cld_scheme=cld_fraction_monan
config_radt_lw_scheme=rrtmg_lw
config_radt_sw_scheme=rrtmg_sw
config_pbl_scheme=bl_shocmf
#config_pbl_scheme=bl_mynn
config_mynn_edmf=0

# ## restart:

# config_do_restart = true

config_pcvol=0
config_cporg=${5}
config_gustf=${6}
config_sub3d=${7}
config_restart=false ####################### PAY ATTENTION HERE 
 
  # iop=GFS # PAY ATTENTION HERE 
  iop=ERA5

  if [ $iop == GFS ] ; then 
    
    #idx=29; fdx=29; im=06; chr=00 ; iy=2024
    idx=3; fdx=3; im=07; chr=00 ; iy=2024
    

    daystep=1 # Duvida
  fi

  if [ $iop == ERA5 ] ; then 
    
   idx=24; fdx=24; im=04; chr=00 ; iy=2024

   daystep=1  
 fi

  # link the initialization file to the right iop
  # link the initialization file to the right IOP
  unlink ${DIR_GRID_SPEC}/${mesh}
  /bin/ln -fs ${DIR_GRID_SPEC}/${mesh}_${iop} ${DIR_GRID_SPEC}/${mesh}

  if [ $mesh == x1.10242 ] ; then # 240km <<<<<<<<<<<
    JobName=SHOCMF240
    Nodes=1
    TasksPerNode=1
    # Ntasks = TasksPerNode * Nodes
    Ntasks=1
    Partition=batch
    config_run_duration=30_00:00:00
    output_interval=24:00:00
    config_dt=1800.0
    config_dt_cu=00:30:00
  fi
  if [ $mesh == x1.40962 ] ; then # 120km <<<<<<<<<<<
    #idx=1; fdx=1; im=03; chr=00
    Nodes=1
    TasksPerNode=128
    # Ntasks = TasksPerNode * Nodes
    Ntasks=128
    Partition=batch
    config_run_duration=2_00:00:00
    output_interval=1:00:00
    config_dt=900.0
    config_dt_cu=00:30:00
  fi
  if [ $mesh == x1.163842 ] ; then # 60km <<<<<<<<<<<
    Nodes=4
    TasksPerNode=128
    # Ntasks = TasksPerNode * Nodes
    Ntasks=512
  
    Partition=PESQ1
    # Partition=batch
    JobName=B_60km
    config_run_duration=10_00:00:00   #PAY ATTENTION HERE
    output_interval=06:00:00
    config_dt=200.0
    config_dt_cu=00:15:00
    time_fila=12:00:00

  fi

  if [ $mesh == x1.256002 ] ; then #  48km <<<<<<<<<<<
    Nodes=2
    TasksPerNode=128
    # Ntasks = TasksPerNode * Nodes
    Ntasks=256
    Partition=batch
    config_run_duration=2_00:00:00
    output_interval=1:00:00
    config_dt=300.0
    config_dt_cu=00:15:00
  fi

  if [ $mesh == x1.655362 ] ; then #  30km <<<<<<<<<<<
    Nodes=4
    TasksPerNode=128
    # Ntasks = TasksPerNode * Nodes
    Ntasks=512
  
    Partition=PESQ1
    # Partition=batch
    JobName=LT1HD050
    config_run_duration=10_00:00:00   #PAY ATTENTION HERE
    output_interval=06:00:00
    config_dt=200.0
    config_dt_cu=00:15:00
    time_fila=12:00:00

  fi
  if [ $mesh == x1.1024002 ] ; then #  24km <<<<<<<<<<<
    #idx=13;fdx=13;im=04
    TasksPerNode=64
    Nodes=16
    # Ntasks = TasksPerNode * Nodes
    Ntasks=1024  
    Partition=batch
    config_run_duration=2_00:00:00
    output_interval=1:00:00
    config_dt=150.0
    config_dt_cu=00:15:00

  fi
  if [ $mesh == x1.2621442 ] ; then #  15km <<<<<<<<<<<
   # idx=26;fdx=26;im=05;iy=2024
    TasksPerNode=80
    Nodes=6
    Ntasks=480
    # Ntasks = TasksPerNode * Nodes  
    Partition=batch
    # Partition=big
    JobName=B_15km
    config_run_duration=06_00:00:00
    output_interval=06:00:00
    config_dt=100.0
    config_dt_cu=00:15:00
    RequestedTime=24:00:00
  fi
  #-- for long runs 
  if [ $monthly_run == 1 ]; then 
     config_run_duration=30_00:00:00
     output_interval=1_00:00:00
  fi 


  #-- for fases 1 or 3
  if [ $fase != 2 ] ; then 
    Nodes=1
    TasksPerNode=128 
    # Ntasks = TasksPerNode * Nodes
    Ntasks=128
    if [ $mesh == x1.10242   ] ; then          TasksPerNode=64 ; Ntasks=64;  fi
    if [ $mesh == x1.40962   ] ; then          TasksPerNode=64 ; Ntasks=64;  fi
    if [ $mesh == x1.2621442 ] ; then Nodes=4; TasksPerNode=64 ; Ntasks=256; fi
    if [ $mesh == x1.655362  ] ; then Nodes=4; TasksPerNode=64 ; Ntasks=256; fi
  fi
#-------------------

  zero=0
  run=1
  nruns=1

  #idx=${3}
  #fdx=${4}
  ini_day=${iy}-${im}-${idx}_00:00:00
  end_day=${iy}-${im}-${fdx}_00:00:00
  echo '===================================================='
  echo 'running MONAN with exec=' ${exec} ${mesh}
  echo 'days=' ${ini_day} ' - ' ${end_day}
  echo 'doing fase:' $fase
  echo '===================================================='

while [ $run -le $nruns ]
do
#======================#======================#======================#======================
 
 new=${exp}.bash
 expname=${exp}_${iop}_${mesh}
 dir_expname=${expname}
 if [ $fase -lt 3 ] ; then
  if [ ! -d ${DIR0}/OLD_RUNS ] ; then mkdir ${DIR0}/OLD_RUNS ; fi
  mkdir  ${DIR0}/${dir_expname}
  unlink ${DIR0}/lastest_run
  /bin/ln -fs ${DIR0}/${dir_expname} ${DIR0}/lastest_run
  #rm -rf ${DIR1}/${dir_expname}
  /bin/cp run.bash  ${expname}.bash
 fi
 cd ${dir_expname} 
 /bin/cp ../load_monan_app_modules.sh .
 
 if [ $fase -lt 5 ] ; then  

 iday=$idx
 while [ $iday -le $fdx ]
 do
   
   cday=$iday
   if [ $cday -lt 10 ] ; then 
     cday=$zero$iday
   fi
   stard=${iy}-${im}-${cday}_${chr}:00:00         # PAY ATTENTION HERE
   echo 'doing ' ${stard} 
   initfiledate=${iy}${im}${cday}${chr}.init.nc
   start_data_dir=${iy}${im}${cday}${chr}
   mkdir ${start_data_dir} 
   cd ${start_data_dir}

   if [ $fase -lt 3 ] ; then
   mkdir hist diag
   /bin/cp ../../ln.sh . 
   ./ln.sh
   /bin/cp ../../nc.clean.sh . 
   /bin/cp ../../clean.sh .
   #----------------------#----------------------
   # GF namelist
   if [ $change_gf_namelist == 1 ] ; then
         /bin/cp ../../GF_ConvPar_nml.bash .
         ./GF_ConvPar_nml.bash 

         sed -i  's/convection_tracer = 0,/convection_tracer = 1,/g' GF_ConvPar_nml
         sed -i  's/use_rhu_ctrl_entr = 0,/use_rhu_ctrl_entr = 0,/g' GF_ConvPar_nml
         sed -i  's/cum_entr_rate     = 1.0e-4, 1.e-3, 2.e-3,/cum_entr_rate   = 8.0e-4, 4.e-3, 1.e-3,/g' GF_ConvPar_nml
         sed -i  's/dicycle           = 0,/dicycle           = 0,/g' GF_ConvPar_nml
         sed -i  's/cum_min_cloud_depth = 200.,100.,150.,/cum_min_cloud_depth = 0.,0.,0.,/g' GF_ConvPar_nml
         sed -i  's/add_coldpool_trig = 0,/add_coldpool_trig = 2,/g' GF_ConvPar_nml 
         sed -i  's/add_coldpool_clos = 0,/add_coldpool_clos = 2,/g' GF_ConvPar_nml 
         sed -i  's/use_pass_cloudvol = 0,/use_pass_cloudvol = 2,/g' GF_ConvPar_nml
         sed -i  's/use_memory        = 0,/use_memory        = 22,/g' GF_ConvPar_nml
         sed -i  's/use_lcl_ctrl_entr = 0,/use_lcl_ctrl_entr = 2,/g' GF_ConvPar_nml

        # sed -i  's/c0_deep           = 1.0e-3,/c0_deep      = 2.0e-3,/g' GF_ConvPar_nml
        # sed -i  's/output_sound      = 0,/output_sound      = 3,/g' GF_ConvPar_nml
        # sed -i  's/closure_choice    = 10,10,3,/closure_choice = 10,1,3,/g' GF_ConvPar_nml
        # sed -i  's/icumulus_gf       = 1,1,0,/icumulus_gf   = 1,1,1,/g' GF_ConvPar_nml 
        # sed -i  's/downdraft         = 1,/downdraft         = 0,/g' GF_ConvPar_nml
        # sed -i  's/use_rebcb         = 1,/use_rebcb         = 0,/g' GF_ConvPar_nml
  fi
  if [ $change_gf_namelist == 0 ] ; then
          if [ -z "GF_ConvPar_nml_${2}.bash" ] ;then
           echo GF_ConvPar_nml_${2}.bash 'does not exist, exiting'
           exit
          fi
          /bin/cp ../../GF_ConvPar_nml_${2}.bash .
          ./GF_ConvPar_nml_${2}.bash .
  fi
  #----------------------#----------------------
 
   # a) streams.init_atm
   /bin/cp ${DIR0}/mesh/${generic_namelists}/streams.init_atmosphere .
   sed -i  's/x1.10242/'${mesh}'/g' streams.init_atmosphere
   sed -i  's/init.nc/'${initfiledate}'/g' streams.init_atmosphere

   # b) namelist.init_atm
   /bin/cp ${DIR0}/mesh/${generic_namelists}/namelist.init_atmosphere .
   sed -i  's/x1.10242/'${mesh}'/g'   namelist.init_atmosphere
   sed -i  's/ini_day/'${stard}'/g' namelist.init_atmosphere
   sed -i  's/end_day/'${stard}'/g' namelist.init_atmosphere
   sed -i  's/CFSR/'${iop}'/g' namelist.init_atmosphere
 
   # c) streams.atm
   /bin/cp ${DIR0}/mesh/${generic_namelists}/streams.atmosphere .
   sed -i  's/x1.10242/'${mesh}'/g' streams.atmosphere
   sed -i  's/init.nc/'${initfiledate}'/g' streams.atmosphere
   sed -i  's/3:00:00/'${output_interval}'/g' streams.atmosphere

   # c) namelist.atm
   /bin/cp ${DIR0}/mesh/${generic_namelists}/namelist.atmosphere .
   sed -i  's/YYYYY/'${stard}'/g' namelist.atmosphere
   sed -i  's/ZZZZZ/'${config_run_duration}'/g' namelist.atmosphere
   sed -i  's/XXXX.X/'${config_dt}'/g' namelist.atmosphere
   sed -i  's/xN.NNNNN/'${mesh}'/g' namelist.atmosphere
   sed -i  's/SSSSS/'${config_physics_suite}'/g' namelist.atmosphere
   sed -i  's/CCCCC/'${config_dt_cu}'/g' namelist.atmosphere
   sed -i  's/CLOUDFRAC/'${config_radt_cld_scheme}'/g' namelist.atmosphere
   sed -i  's/_LW_/'${config_radt_lw_scheme}'/g' namelist.atmosphere
   sed -i  's/_SW_/'${config_radt_sw_scheme}'/g' namelist.atmosphere
   sed -i  's/_PBL_/'${config_pbl_scheme}'/g' namelist.atmosphere
   sed -i  's/PCVOL/'${config_pcvol}'/g' namelist.atmosphere
   sed -i  's/CPORG/'${config_cporg}'/g' namelist.atmosphere
   sed -i  's/GUSTF/'${config_gustf}'/g' namelist.atmosphere
   sed -i  's/SUB3D/'${config_sub3d}'/g' namelist.atmosphere
   sed -i  's/ZERO/'${config_mynn_edmf}'/g' namelist.atmosphere
   sed -i  's/restart_false/'${config_restart}'/g' namelist.atmosphere

   if [ $mesh == x1.1024002 ] ; then
       sed -i  's/config_len_disp = 0.0/config_len_disp = 24000.0/g' namelist.atmosphere
   fi
   # d) copy the output namelists 
   /bin/cp ${DIR0}/mesh/${generic_namelists}/stream_list.atmosphere.surface  . 
   /bin/cp ${DIR0}/mesh/${generic_namelists}/stream_list.atmosphere.output  . 
   /bin/cp ${DIR0}/mesh/${generic_namelists}/stream_list.atmosphere.diagnostics  . 
   
   ### only for diurnal cycle
   /bin/cp ${DIR0}/mesh/${generic_namelists}/stream_list.atmosphere.diagnostics_dc  . 

   fi # fase < 3

   #---------------------------------------------------------------------------------------
   # d) run model /initialization
   if [ $fase -lt 3 ] ; then 
      twentyfour=24
      ndays=${config_run_duration%_*}
      outint=${output_interval%:00:00*}
      step=$(( $twentyfour /  $outint ))
      noutputs=$(( 1 + $ndays * $step ))
      #mpirun -np 1 ${exec} > model_${exp}.out &
      #template_submit.bash ${Nodes} ${Partition} ${TasksPerNode} ${Ntasks} ${exec}
      #chmod +x submit.bash
      #/usr/bin/sbatch submit.bash 
      submit=submit_${expname}_${start_data_dir}.bash
#---------------------------------------------------------------------------------------

cat << Eof > ${submit}
#!/bin/bash
#SBATCH --job-name=${JobName}
#SBATCH --nodes=${Nodes}
#SBATCH --partition=${Partition}
#SBATCH --tasks-per-node=${TasksPerNode}
#SBATCH --ntasks=${Ntasks}
#SBATCH --time=${RequestedTime}  
#SBATCH --output=/home/guilherme.farache/logs/%j.log
#SBATCH --exclusive

      cd \$SLURM_SUBMIT_DIR
      echo \$SLURM_SUBMIT_DIR

      module purge
      module load ohpc
      module unload openmpi4
      module load phdf5
      module load netcdf
      module load netcdf-fortran
      module load mpich-4.0.2-gcc-9.4.0-gpof2pv
      module load hwloc
      module load imagemagick-7.0.8-7-gcc-11.2.0-46pk2go 

      export OMP_NUM_THREADS=1
      export OMPI_MCA_btl_openib_allow_ib=1
      export OMPI_MCA_btl_openib_if_include="mlx5_0:1"
      export PMIX_MCA_gds=hash
      export NETCDF=/mnt/beegfs/monan/libs/netcdf
      export PNETCDF=/mnt/beegfs/monan/libs/PnetCDF

      MPI_PARAMS="-iface ib0 -bind-to core -map-by core"
      export MKL_NUM_THREADS=1
      export I_MPI_DEBUG=5
      export MKL_DEBUG_CPU_TYPE=5
      export I_MPI_ADJUST_BCAST=12 ## NUMA aware SHM-Based (AVX512)
      # PIO is not necessary for version 8.* If PIO is empty, MPAS Will use SMIOL
      export PIO=
      export PMIX_MCA_gds=hash
      ulimit -c unlimited
      ulimit -v unlimited
      ulimit -s unlimited
     
      ###mpirun -genvall ${exec} &> ${exec}.log 
      mpirun -env UCX_NET_DEVICES=mlx5_0:1 -genvall ${exec} &> ${exec}.log
      
      #----------------------
      export fase=$fase
      export doconvert=1 ; export local_time=1 ; export number_days=${ndays}; export noutputs=${noutputs}; export step=${step}
      cd ${DIR0}/${dir_expname}/${start_data_dir}
      /bin/cp ../../process_monan_fase1.bash ${new}
      ${new} ${expname} ${mesh} ${start_data_dir} ${idx} ${fdx} ${im} ${iy} ${cday} ${chr} ${initfiledate} >> ${exec}.log &
      ### dc 
      #/bin/cp ../../process_monan_fase1_dc.bash ${new}
      #${new} ${expname} ${mesh} ${start_data_dir} ${idx} ${fdx} ${im} ${iy} ${cday} ${chr} ${initfiledate} >> ${exec}.log


Eof

#--------------------------------------------------------------------------------------- 
      chmod +x ${submit}
      /usr/bin/sbatch ${submit}
      echo 'submit' $start_data_dir $submit
      #sleep 10
   fi
   # e) pos-proc 1
   if [ $fase == 3 ] || [ $fase == 4 ] ; then
      export fase=$fase
      twentyfour=24
      ndays=${config_run_duration%_*}
      outint=${output_interval%:00:00*}
      step=$(( $twentyfour /  $outint ))
      noutputs=$(( 1 + $ndays * $step ))
      
      export doconvert=1 ; export local_time=1 ; export number_days=${ndays}; export noutputs=${noutputs}; export step=${step}
      cd ${DIR0}/${dir_expname}/${start_data_dir}
      /bin/cp ../../process_monan_fase1.bash ${new}
      ${new} ${expname} ${mesh} ${start_data_dir} ${idx} ${fdx} ${im} ${iy} ${cday} ${chr} ${initfiledate} > process.log
      # dc 
      /bin/cp ../../process_monan_fase1_dc.bash ${new}
      ${new} ${expname} ${mesh} ${start_data_dir} ${idx} ${fdx} ${im} ${iy} ${cday} ${chr} ${initfiledate} >> ${exec}.log
  
   fi 

   #--------------------------------
   # back to the model directory 
   cd ${DIR0}/${dir_expname}
   #--------------------------------
  # next day
  iday=`expr $iday + $daystep`
done
fi # fase < 5

# f) pos-proc 2
if [ $fase == 5 ] || [ $fase == 6 ] ; then
#  [ $x == "Y" ] || [ $x == "y"]
      
      cd ${DIR0}/${dir_expname}
      /bin/cp ../process_monan_fase2.bash .
      process_monan_fase2.bash ${expname} ${fase} ${mesh} ${idx} ${fdx} ${im} ${iy} ${chr} ${start_data_dir} 
fi 

run=`expr $run + 1`

done
exit





