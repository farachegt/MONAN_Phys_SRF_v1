#!/bin/bash
DIR0=/mnt/beegfs/saulo.freitas/runs/MONAN/model
DIR1=/mnt/beegfs/saulo.freitas/runs/MONAN/posprocess
DIR_GRID_SPEC=/mnt/beegfs/saulo.freitas/runs/MONAN/model/datain
#DIR_CONVERT_MPAS=/home/saulo.freitas/models/MONAN/convert_mpas_pkubota/convert_mpas_v0.1.0_egeon.gnu940/orig
DIR_CONVERT_MPAS=/home/saulo.freitas/models/MONAN/convert_mpas_pkubota/convert_mpas_v0.1.0_egeon.gnu940/vizmp
DIR_EVAL=/mnt/beegfs/saulo.freitas/runs/MONAN/eval_data
DIR_GRADS=/home/saulo.freitas/bin/opengrads-2.2.1.oga.1/Contents
local_shells=/mnt/beegfs/saulo.freitas/runs/scripts_shells
DIR_imagemagick=/home/saulo.freitas/bin/imagemagick/bin
#-----------------------------------------------

if [ $fase -lt 2 ] ; then exit; fi

../load_monan_app_modules.sh
module load phdf5
module load netcdf
module load netcdf-fortran
export NETCDF=/mnt/beegfs/monan/libs/netcdf
export PNETCDF=/mnt/beegfs/monan/libs/PnetCDF
module load imagemagick-7.0.8-7-gcc-11.2.0-46pk2go 
##module load cdo-2.0.4-gcc-9.4.0-bjulvnd
maxpostpernode=31    # <------ qtde max de convert_mpas por no!

#--------------------- input data
expname=${1}
mesh=${2}
start_data_dir=${3}
idx=${4}
fdx=${5}
im=${6}
iy=${7}
cday=${8}
chr=${9}
initfiledate=${10}
#domain=${11}

echo $DIR0 $expname  ${initfiledate} ${start_data_dir}
#--------------------- diurnal cycle processing (per fcst day).
cd ${DIR0}/${expname}/${start_data_dir}/diag
/bin/cp ${DIR1}/pos_config_${mesh}/* .
/bin/ln -fs target_domain_${domain} target_domain
/bin/ln -fs ${DIR_GRID_SPEC}/${domain}/${mesh}/${mesh}.${initfiledate} file_grid_spc.init.nc
/bin/ln -fs ${DIR_CONVERT_MPAS}/convert_mpas .
/bin/cp ${local_shells}/convert_mpas_grads.bash .
 #--- to convert to local time
 /bin/cp ${local_shells}/monan_conv.gs conv.gs
 /bin/cp ${local_shells}/${mesh}.ctl x.ctl
 /bin/cp ${local_shells}/lt.inp .
 /bin/cp ${local_shells}/lt.x .
 /bin/cp ${local_shells}/LT_${mesh}.ctl LT_x.ctl

#---------------------
zero=0
mkdir fct1 fct2 fct3 fct4 fct5 backup backup2
rm -f latlon.nc
var2d=ctt,u_gustfront,v_gustfront,uzonal_850hPa,umeridional_850hPa,wlpool,vcpool,v2d1,rmfxshcu,rtopdpcu,\
rtopmdcu,rtopshcu,rbotdpcu,rmfxmdcu,rainc,rainnc,cuprec,gsw,glw,hpbl,cldfrac_tot_UPP,t2m,lh,hfx,sst,precipw,olrtoa,cape,cin,precipci

var2d_aval=ctt,gsw,lh,hfx,t2m,rmfxdpcu,rmfxshcu,cldfrac_tot_UPP,rainnc,rainc,rmfxmdcu,\
relhum_250hPa,relhum_500hPa,relhum_850hPa,relhum_950hPa,relhum_150hPa,relhum_300hPa,relhum_700hPa,\
temperature_250hPa,temperature_500hPa,temperature_850hPa,temperature_950hPa,temperature_150hPa,temperature_300hPa,temperature_700hPa,\
qv_250hPa,qv_500hPa,qv_850hPa,qv_950hPa,qv_150hPa,qv_300hPa,qv_700hPa,\
zgeo_250hPa,zgeo_500hPa,zgeo_850hPa,zgeo_950hPa,zgeo_150hPa,zgeo_300hPa,zgeo_700hPa,\
rthcuten_200hPa,rthcuten_250hPa,rthcuten_300hPa,rthcuten_500hPa,rthcuten_700hPa,rthcuten_850hPa,rthcuten_950hPa,\
rqvcuten_200hPa,rqvcuten_250hPa,rqvcuten_300hPa,rqvcuten_500hPa,rqvcuten_700hPa,rqvcuten_850hPa,rqvcuten_950hPa,\
rqvblten_700hPa,rqvblten_825hPa,rqvblten_875hPa,rqvblten_900hPa,rqvblten_700hPa,rqvblten_850hPa,rqvblten_950hPa,\
rthblten_700hPa,rthblten_825hPa,rthblten_875hPa,rthblten_900hPa,rthblten_700hPa,rthblten_850hPa,rthblten_950hPa,\
rthcuten_400hPa,rthcuten_600hPa,rthcuten_825hPa,\
rqvcuten_400hPa,rqvcuten_600hPa,rqvcuten_825hPa

var3d=pressure,cldfrac,qv,qc,qi,rthcuten,rqvcuten,rqccuten,rqicuten,rdnmfxcu,rupmfxcu,v3d1
extf=${expname}-${iy}-${im}-${cday}

diag_list=()
diag_dc_list=()
rm -f latlon.nc

#------------------#------------------
if [ ${doconvert} -ge 1 ] ; then 

  if [[ -e backup ]]; then
    /bin/mv backup/diag_dc.* .
  fi

  #--- get the list of files to be converted 
  x=0
  #--- loop 1 
  for f in diag_dc.*
  do
    x=`expr $x + 1`
    diag_dc_list[${x}]=${f}
    echo 'f=' ${diag_dc_list[$x]}

  done  #--- end of the loop 1
  i=1
  f1=all_${diag_dc_list[$i]}
  f2=only_precip_${diag_dc_list[$i]}
  f3=only_2d_aval_${diag_dc_list[$i]}

  if [[ -e ${diag_dc_list[1]} ]] && [[ ${doconvert} == 1 ]]; then 

    #--- get the convert the files 
    node=1
    inicio=1
    nfiles=${x}
    fim=$((maxpostpernode <= nfiles ? maxpostpernode : nfiles))
  
    while [ ${inicio} -le ${nfiles} ]
      do
       echo '===>'  ${inicio} ${nfiles} 
      for i in $(seq  ${inicio} ${fim})
       do
        echo 'doing ==>' $i ${diag_dc_list[$i]}
        ./convert_mpas file_grid_spc.init.nc ${diag_dc_list[$i]} ll_${diag_dc_list[$i]} >& conv_diag_dc_${i}.out &
       done
      inicio=$((fim + 1))
      temp=$((fim + maxpostpernode))
      fim=$(( temp < nfiles ? temp : nfiles ))
      node=$((node+1))
      # necessario aguardar as rodadas em background
      wait
    done

    #--- 
    #
    #-- get the hourly precip
    inicio=1
    /bin/cp  ll_${diag_dc_list[1]} previous_precip.nc
    for i in $(seq  ${inicio} ${fim})
      do
         echo 'precip between: '$i ll_${diag_dc_list[$i]} previous_precip.nc
         cdo -L sub -apply,-selname,rainc,rainnc [ ll_${diag_dc_list[$i]} previous_precip.nc ] ll_precip_${i}.nc
         #-- save for next step
         /bin/cp ll_${diag_dc_list[$i]} previous_precip.nc
         cdo cat ll_precip_${i}.nc ${f2}
         /bin/mv ll_precip_${i}.nc  backup2
    done
    #---
    #-- for troposphere evaluation
    for i in $(seq  ${inicio} ${fim})
      do          
        #-- put all files into only one  ${f1} and ${f2} 
        cdo selname,${var2d_aval} ll_${diag_dc_list[$i]} ${i}_2d_aval_ll.nc
        cdo cat ${i}_2d_aval_ll.nc  ${f3}
        cdo cat ll_${diag_dc_list[$i]}  ${f1}

        #-- move processed files to backup2
        /bin/mv ll_${diag_dc_list[$i]} ${i}_2d_aval_ll.nc   backup2
        #-- move original files to backup
        /bin/mv ${diag_dc_list[$i]} backup
    done
  fi  
  inicio=1
  for i in $(seq  ${inicio} ${fim})
    do       
       if [[ -e ${diag_dc_list[$i]} ]] ; then  
           /bin/mv ${diag_dc_list[$i]} backup
       fi 
  done 
  mkdir conv.out
  /bin/mv con*.out conv.out
#------------------
#------------------
  #--- getting the mean daily to compare with era-5 or other reference data.
  /usr/bin/unlink monan_zqtr_data.nc

  if [ ${number_of_days} == 30  ] ; then 
      cdo monmean ${f3} timemean_${f3}
      export number_of_days=1
  else
      cdo daymean ${f3} timemean_${f3}
  fi
  #--- getting the initial state for removing bias from the fct
  cdo seltimestep,1  ${f3} ci_${f3}
  cdo -r settaxis,${iy}-${im}-${cday},${chr}:00:00,1day ci_${f3} x_ci_${f3}
  /bin/mv  x_ci_${f3} ci_${f3}
  /bin/ln -fs ci_${f3} ci_monan_zqtr_data.nc
  
  #--- part 1
   cdo -r settaxis,${iy}-${im}-${cday},${chr}:00:00,1day timemean_${f3} xtimemean_${f3}
   /bin/mv xtimemean_${f3} timemean_${f3}
   /bin/ln -fs timemean_${f3} monan_zqtr_data.nc
   
   #--- part 2
   #--- converting era5 to MONAN grid configuration (casoRS), including data used for MONAN's IC
    echo 'evaluating case=' ${case}
   if [ ${case} == casoRS  ] ; then export era5_data=era5_ave_00Z224_23Z30_042024.nc; export ci_era5_data=ci_era5_00Z24042024.nc; fi
   if [ ${case} == aug2023 ] ; then export era5_data=era5_ave_00Z225_23Z31_082023.nc; export ci_era5_data=ci_era5_00Z25082023.nc; fi
   if [ ${case} == jan2024 ] ; then export era5_data=mon_mean_jan2024_era5.nc; fi
   if [ ${case} == jul2024 ] ; then export era5_data=mon_mean_jul2024_era5.nc; fi

   /bin/cp ${DIR_EVAL}/ERA5/${case}/${era5_data} .
   cdo -remapycon,timemean_${f3} ${era5_data} era5_zqtr_data.nc
   
   /bin/cp ${DIR_EVAL}/ERA5/${case}/${ci_era5_data} .
   cdo -remapycon,timemean_${f3} ${ci_era5_data} ci_era5_zqtr_data.nc
   
   #--- precip 
   /bin/rm previous_precip.nc
   # fixing the time interval
   #cdo -r settaxis,2000-01-01,12:00:00,1h x.nc y.nc
   #cdo -r settaxis,${iy}-${im}-${cday},${chr}:00:00,1h ${f1} x_${f1}
   #cdo -r settaxis,${iy}-${im}-${cday},${chr}:00:00,1h ${f2} x_${f2}
   #/bin/mv x_${f1} ${f1}
   #/bin/mv x_${f2} ${f2}

   #-- getting the mean diurnal cycle,  fixing the time interval
   cdo dhourmean ${f1} dc_${f1}
   cdo dhourmean ${f2} dc_${f2}
   cdo -r settaxis,${iy}-${im}-${cday},${chr}:00:00,1h dc_${f1} x_dc_${f1}
   cdo -r settaxis,${iy}-${im}-${cday},${chr}:00:00,1h dc_${f2} x_dc_${f2}
   /bin/mv x_dc_${f1} dc_${f1}
   /bin/mv x_dc_${f2} dc_${f2}
  
   #--- link for pngs 
   /bin/ln -fs dc_${f2} precip_dc.nc

fi # doconvert loop 
###############################################################

#--- pngs 
#----------------------------------------------------------------------------------------------------------------------
#
#--- evaluation with ERA-5
#--- create grads script/generate pngs
/bin/rm *.aval.png *.aval-?.png *_ave*day_*.png *aval_mean-?.png *aval_mean.png
mkdir pngs

#------------------
ff=${number_of_days}
cat << Eof3 > tmp3.gs
    ff='${ff}'
    title='${expname}'
    biasrm='${biasrm}'
Eof3
cat ./tmp3.gs ${local_shells}/mean_aval_monan.gs > mean_p3.gs 
${DIR_GRADS}/grads -blc 'run mean_p3.gs' >mean_p3.out 

#--- montages
ivar=1 ; nvar=9
while [ $ivar -le $nvar ]
do
if [ $ivar == 1 ] ; then var=TEMP ;fi
if [ $ivar == 2 ] ; then var=QV ;fi
if [ $ivar == 3 ] ; then var=RH ;fi
if [ $ivar == 4 ] ; then var=Z ;fi
if [ $ivar == 5 ] ; then var=TTCU ;fi
if [ $ivar == 6 ] ; then var=TQCU ;fi
if [ $ivar == 7 ] ; then var=TTBL ;fi
if [ $ivar == 8 ] ; then var=TQBL ;fi
if [ $ivar == 9 ] ; then var=X ;fi

ireg=1; nreg=2
while [ $ireg -le $nreg ]
do
 if [ $ireg == 1 ] ; then reg=SAM ;fi
 if [ $ireg == 2 ] ; then reg=GLB ;fi
 if [ $ireg == 3 ] ; then reg=RSU ;fi

 echo $reg $var
 fig1=${reg}_${var}
 export fcst_files=`/bin/ls -1tr ${fig1}*.png`
 montage -tile 3x3 -geometry +1+1 -background white $fcst_files  ${fig1}.aval_mean.png
 /bin/mv $fcst_files pngs

ireg=`expr $ireg + 1`
done
ivar=`expr $ivar + 1`
done
 
/usr/bin/tar -cvf ${expname}_${start_data_dir}.aval_mean.tar *.aval_mean.png *_ave*day*.png
/usr/bin/mv ${expname}_${start_data_dir}.aval_mean.tar /home/saulo.freitas
/bin/mv *.png pngs

#------------------
exit
#------------------

cat << Eof3 > tmp4.gs
    ff='${ff}'
    title='${expname}'
Eof3
cat ./tmp4.gs ${local_shells}/aval_monan.gs > p3.gs 
${DIR_GRADS}/grads -blc 'run p3.gs' >p3.out 
#--- montages
ivar=0 ; nvar=22
while [ $ivar -le $nvar ]
do
if [ $ivar == 0 ] ; then var=TEMP2M ;fi
if [ $ivar == 1 ] ; then var=TEMP950hPa ;fi
if [ $ivar == 2 ] ; then var=TEMP850hPa ;fi
if [ $ivar == 3 ] ; then var=TEMP500hPa ;fi
if [ $ivar == 4 ] ; then var=TEMP250hPa ;fi
if [ $ivar == 5 ] ; then var=QV950hPa ;fi
if [ $ivar == 6 ] ; then var=QV850hPa ;fi
if [ $ivar == 7 ] ; then var=QV500hPa ;fi
if [ $ivar == 8 ] ; then var=QV250hPa ;fi
if [ $ivar == 9 ] ; then var=RH950hPa ;fi
if [ $ivar == 10 ] ; then var=RH850hPa ;fi
if [ $ivar == 11 ] ; then var=RH500hPa ;fi
if [ $ivar == 12 ] ; then var=RH250hPa ;fi
if [ $ivar == 13 ] ; then var=Z950hPa ;fi
if [ $ivar == 14 ] ; then var=Z850hPa ;fi
if [ $ivar == 15 ] ; then var=Z500hPa ;fi
if [ $ivar == 16 ] ; then var=Z250hPa ;fi
if [ $ivar == 17 ] ; then var=DPMF ;fi
if [ $ivar == 18 ] ; then var=SHMF ;fi
if [ $ivar == 19 ] ; then var=CFRA ;fi
if [ $ivar == 20 ] ; then var=GSW ;fi
if [ $ivar == 21 ] ; then var=LE ;fi
if [ $ivar == 22 ] ; then var=H ;fi

ireg=1; nreg=2
while [ $ireg -le $nreg ]
do
 if [ $ireg == 1 ] ; then reg=SAM ;fi
 if [ $ireg == 2 ] ; then reg=GLB ;fi
 if [ $ireg == 3 ] ; then reg=RSU ;fi

 echo $reg $var
 fig1=${reg}_${var}
 export fcst_files=`/bin/ls -1tr ${fig1}*.png`
 if [[ ${ff} -lt 7 ]] ; then 
    #${DIR_imagemagick}/montage -tile 3x2 -geometry +1+1 -background white $fcst_files  ${fig1}.aval.png
     montage -tile 3x2 -geometry +1+1 -background white $fcst_files  ${fig1}.aval.png

 fi
 if [[ ${ff} -gt 6 ]] ; then 
   #${DIR_imagemagick}/montage -tile 3x3 -geometry +1+1 -background white $fcst_files  ${fig1}.aval.png
                       montage -tile 3x3 -geometry +1+1 -background white $fcst_files  ${fig1}.aval.png
 fi

 ireg=`expr $ireg + 1`
done
ivar=`expr $ivar + 1`
done
/usr/bin/tar -cvf ${expname}_${start_data_dir}.aval.tar *.aval.png
/usr/bin/mv ${expname}_${start_data_dir}.aval.tar /home/saulo.freitas
mkdir pngs; /bin/mv *.png pngs


#----------------------------------------------------------------------------------------------------------------------
exit
#----------------------------------------------------------------------------------------------------------------------



#--- precip evaluation
#--- link monan data
export monandata=`/bin/ls -1 dc_only_precip*`
/bin/ln -fs  ${monandata} precip_dc.nc

#--- UTC TIME
if [ ${local_time} == 0 ] ; then
   title2=UTC
   pref=
fi
#--- LOCAL TIME
if [ ${local_time} == 1 ] ; then
   title2=LT
   pref=LT_
   #-- transform MONAN data from UTC to LOCAL TIME
   /bin/ln -fs  ${monandata} x.nc
   ${DIR_GRADS}/grads -blc 'run conv.gs' >conv.out 
   ./lt.x > lt1.out
   cdo -f nc import_binary LT_x.ctl LT_x.nc
   /bin/mv LT_x.nc LT_${monandata}
   /bin/ln -fs LT_${monandata} precip_dc.nc
fi
#--- data time period
# MAR 2024
# 30 days
#timeperiod=202403_00Z01_24Z31
# 10 days
#timeperiod=202403_00Z01_24Z10
# JAN 2024
timeperiod=202401_00Z07_24Z16
#title2=LT; pref=LT_

#--- link land_mask
f=land_mask.nc
if ! [ -e ${f} ]; then 
#/bin/ln -fs ${DIR_EVAL}/land_mask/${pref}${mesh}_land_mask.nc land_mask.nc
#--- creating the land_mask for this run
/bin/rm -f land_mask.nc
/bin/cp ${DIR_EVAL}/land_mask/era5-land-mask.nc .
cdo -remapycon,precip_dc.nc era5-land-mask.nc land_mask.nc
fi

#--- get ERA5 data
f=monan_era5.nc4
if ! [ -e ${f} ]; then 
/bin/cp ${DIR_EVAL}/ERA5/${pref}ERA5.${timeperiod}_ave_dc1h.nc era5.nc4
cdo -remapycon,precip_dc.nc era5.nc4 monan_era5.nc4
fi

#--- get GPM data
f=monan_gpm.nc4
if ! [ -e ${f} ]; then 
/bin/cp ${DIR_EVAL}/GPM/${pref}3B-HHR-L.MS.MRG.3IMERG.${timeperiod}_ave_dc1h.nc4 gpm.nc4
cdo -remapycon,precip_dc.nc gpm.nc4 monan_gpm.nc4
fi

#--- get GPM clim 2021/2022 JAN
f=monan_gpm_clim.nc4
if ! [ -e ${f} ]; then 
/bin/cp ${DIR_EVAL}/GPM/LT_3B-HHR-L.MS.MRG.3IMERG.2022and2021_01_00Z07_24Z16_ave_dc1h.nc4 gpm_clim.nc4
#/bin/cp ${DIR_EVAL}/GPM/LT_3B-HHR-L.MS.MRG.3IMERG.2022and2021_01_00Z01_24Z30_ave_dc1h.nc4 gpm_clim.nc4
cdo -remapycon,precip_dc.nc gpm_clim.nc4 monan_gpm_clim.nc4
fi


#--- link GFS data
#/bin/ln -fs ${DIR_EVAL}/GFS/${pref}${mesh}_GFS_${timeperiod}.nc monan_gfs.nc4

#--- create grads script/generate pngs
cat << Eof3 > tmp.gs
    title='${expname}'
    title2='${title2}'
Eof3
cat ./tmp.gs ${local_shells}/show_dc_monan.gs > p2.gs 

${DIR_GRADS}/grads -blc 'run p2.gs' >p2.out 


#--- montage pngs 
${DIR_imagemagick}/montage -tile 3x3 -border 3 -bordercolor black -geometry +1+1 -background white \
GLOBAL_MONAN.png LAND_MONAN.png OCEANS_MONAN.png \
AMAZON_MONAN.png NEB_MONAN.png SEB_MONAN.png     \
CONUS_MONAN.png AFRICA_MONAN.png SAHEL_MONAN.png tar_${pref}_all1_dc.png

${DIR_imagemagick}/montage -tile 2x3 -border 3  -bordercolor black -geometry +1+1 -background white \
AUST_MONAN.png  INDIA_MONAN.png INDC_MONAN.png  \
EUR_MONAN.png  PACIFIC1_MONAN.png tar_${pref}_all2_dc.png

#---- pngs files part final
mkdir pngs; /bin/mv x_*.png pngs ; /bin/mv *MONAN.png pngs


/usr/bin/tar -cvf ${expname}.tar tar_*.png pngs
#if [ -e "$1" ] ;then
#    echo "file exist"
#    exit
#fi


/usr/bin/mv ${expname}.tar /home/saulo.freitas
mkdir pngs; /bin/mv *.png pngs

mkdir lixo; /bin/mv *.out *.bin lixo


###############################################################
exit 


