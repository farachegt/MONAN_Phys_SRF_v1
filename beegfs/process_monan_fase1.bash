#!/bin/bash
DIR0=/mnt/beegfs/saulo.freitas/runs/MONAN/model
DIR1=/mnt/beegfs/saulo.freitas/runs/MONAN/posprocess
DIR_GRID_SPEC=/mnt/beegfs/saulo.freitas/runs/MONAN/model/datain
DIR_CONVERT_MPAS=/home/saulo.freitas/models/MONAN/convert_mpas_pkubota/convert_mpas_v0.1.0_egeon.gnu940/vizmp
local_shells=/mnt/beegfs/saulo.freitas/runs/scripts_shells
DIR_GRADS=/home/saulo.freitas/bin/opengrads-2.2.1.oga.1/Contents
#DIR_imagemagick=/home/saulo.freitas/bin/imagemagick/bin
DIR_imagemagick=/opt/spack/opt/spack/linux-rhel8-zen2/gcc-11.2.0/imagemagick-7.0.8-7-46pk2gowitva2xtvogtfg2ubug2d6v2w/bin
#-----------------------------------------------

if [ $fase -lt 2 ] ; then exit; fi

../load_monan_app_modules.sh
module load phdf5
module load netcdf
module load netcdf-fortran
export NETCDF=/mnt/beegfs/monan/libs/netcdf
export PNETCDF=/mnt/beegfs/monan/libs/PnetCDF
#smodule load cdo-2.0.4-gcc-9.4.0-bjulvnd
module load imagemagick-7.0.8-7-gcc-11.2.0-46pk2go 

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
/bin/ln -fs ${DIR_GRID_SPEC}/${domain}/${mesh}/${mesh}.${initfiledate} file_grid_spc.init.nc
/bin/ln -fs ${DIR_CONVERT_MPAS}/convert_mpas .
/bin/cp ${local_shells}/convert_mpas_grads.bash .
###############################################################
zero=0
mkdir fct1 fct2 fct3 fct4 fct5 backup 
rm -f latlon.nc
var2d=u_gustfront,v_gustfront,uzonal_850hPa,umeridional_850hPa,wlpool,vcpool,v2d1,rmfxshcu,rtopdpcu,rtopmdcu,rtopshcu,rbotdpcu,rmfxmdcu,rainc,rainnc,cuprec,gsw,glw,hpbl,cldfrac_tot_UPP,t2m,hl,hfx,sst,precipw,olrtoa,cape,cin 
var3d=pressure,cldfrac,qv,qc,qi,rthcuten,rqvcuten,rqccuten,rqicuten,rdnmfxcu,rupmfxcu,cuprec,buoyx,wlpool,rainc,rainnc,v2d1,rtopdpcu,rtopmdcu,rtopshcu,rmfxmdcu,rmfxshcu,rmfxdpcu,precipci,ctt
extf=${expname}-${iy}-${im}-${cday}


if [ $fase -lt 4 ] ; then 

if [ ${doconvert} == 1 ] ; then 

x=0
#--- loop 
for f in diag.*
do

echo 'diag/hist file=' ${f}
./convert_mpas file_grid_spc.init.nc ${f}
/bin/mv latlon.nc ll_${f}
/bin/mv ${f} backup

#cdo selname,u_gustfront,v_gustfront,uzonal_850hPa,umeridional_850hPa,wlpool,vcpool,v2d1,cuprec,rmfxshcu,rtopdpcu,rtopmdcu,rtopshcu,rbotdpcu,rmfxmdcu,rainc,rainnc ll_$f 2d_coldpool_ll_$f
#cdo selname,cuprec,gsw,glw,hpbl,cldfrac_tot_UPP,rainc,rainnc,t2m,hl,hfx,sst,precipw,olrtoa,cape,cin ll_$f 2d_ext_ll_$f
#cdo -L -sellevel,1/42 -selname,pressure,cldfrac,qv,qc,qi,rthcuten,rqvcuten,rqccuten,rqicuten,rdnmfxcu,rupmfxcu ll_$f 3d_ext_ll_$f

cdo selname,${var2d} ll_$f 2d_ext_ll_$f
cdo -L -sellevel,0/42 -selname,${var3d} ll_$f 3d_ext_ll_$f

#cdo select,levrange,25,900,name=T icon_oce.nc sellevel.nc

#### for the diurnal cycle
file=$x
if [[ $x -lt 10 ]] ; then 
file=$zero$x
fi
echo 'x=' $file

file=$x
if [[ $x -lt 10 ]] ; then 
file=$zero$x
fi
echo 'x=' $file

if [[ $x -lt 24 ]] ; then
echo 'moving day1' $x ll_$f 
/bin/mv ll_$f fct1/${file}_ll_$f
fi
#--- copying the 00Z of the next day 
if [[ $x -eq 24 ]] ; then
echo 'copying day1' $x ll_$f 
/bin/cp ll_$f fct1/${file}_ll_$f
fi

if [[ $x -lt 48 ]] && [[ $x -gt 23 ]] ; then 
echo 'moving day2' $x ll_$f 
/bin/mv ll_$f fct2/${file}_ll_$f
fi
#--- copying the 00Z of the next day 
if [[ $x -eq 48 ]] ; then
echo 'copying day2' $x ll_$f 
/bin/cp ll_$f fct2/${file}_ll_$f
fi

if [[ $x -lt 72 ]] && [[ $x -gt 47 ]] ; then 
echo 'moving day3' $x ll_$f 
/bin/mv ll_$f fct3/${file}_ll_$f
fi

#--- copying the 00Z of the next day 
if [[ $x -eq 72 ]] ; then
echo 'copying day3' $x ll_$f 
/bin/cp ll_$f fct3/${file}_ll_$f
fi

if [[ $x -lt 96 ]] && [[ $x -gt 71 ]] ; then 
echo 'moving day4' $x ll_$f 
/bin/mv ll_$f fct4/${file}_ll_$f
fi

#--- copying the 00Z of the next day 
if [[ $x -eq 96 ]] ; then
echo 'copying day4' $x ll_$f 
/bin/cp ll_$f fct4/${file}_ll_$f
fi

if [[ $x -lt 120 ]] && [[ $x -gt 95 ]] ; then 
echo 'moving day5' $x ll_$f 
/bin/mv ll_$f fct5/${file}_ll_$f
fi

#--- copying the 00Z of the next day 
if [[ $x -eq 120 ]] ; then
echo 'copying day5' $x ll_$f 
/bin/cp ll_$f fct5/${file}_ll_$f
fi
### end of the diurnal cycle

x=`expr $x + 1`
done
#--- end of the loop


# put all files into only one
file2d=all_2d_ext_ll_${extf}.nc
if [ -f $file2d ] ; then /bin/rm $file2d; fi
file3d=all_3d_ext_ll_${extf}.nc
if [ -f $file3d ] ; then /bin/rm $file3d; fi
 

export fcst_files=`/bin/ls -1 2d_ext_ll_*.nc`
echo $fcst_files  
cdo cat $fcst_files ${file2d}
/bin/ln -fs  ${file2d} all_2d.nc

export fcst_files=`/bin/ls -1 3d_ext_ll_*.nc`
echo $fcst_files  
cdo cat $fcst_files  ${file3d}
/bin/ln -fs  ${file3d} all_3d.nc

fi # doconvert

#number_days=`expr $x - 1`
#export number_days=${ndays}; export noutputs=${noutputs}; export step=${step}
#number_days=10
#step=2

##fi # fase < 4


 
#---- pngs files part 1
cat << Eof3 > tmp.gs
 ff='${noutputs}'
 step='${step}'
 file1='all_2d'
 title='${expname}'
 title2='${expname} \mean diurnal cycle (mm/day) FCST day 2'
Eof3

#--- mean fields
 cat ./tmp.gs ${local_shells}/mean_precip_monan.gs > mean_p1.gs 
 ${DIR_GRADS}/grads -blc 'run mean_p1.gs' >mean_p1.out 
 export fcst_files=`/bin/ls -1tr *mean*.png`
 montage -tile 3x2 -geometry +1+1 -background white $fcst_files  mean_precips_${expname}_${start_data_dir}.png
 /bin/cp mean_precips_${expname}_${start_data_dir}.png /home/saulo.freitas

#-- daily fields
 cat ./tmp.gs ${local_shells}/precip_monan.gs > p1.gs 
 ${DIR_GRADS}/grads -blc 'run p1.gs' >p1.out 


fi # fase < 4

if [ -d pngs ] ; then
   /bin/mv ./pngs/*.png .
else
   mkdir pngs
fi
/bin/rm *.x*.png

#--- montage pngs 
ivar=1
nvar=7
while [ $ivar -le $nvar ]
do
if [ $ivar == 1 ] ; then var=NCV ;fi
if [ $ivar == 2 ] ; then var=CNV ;fi
if [ $ivar == 3 ] ; then var=TOT ;fi
if [ $ivar == 4 ] ; then var=CLDFRTOT ;fi
if [ $ivar == 5 ] ; then var=T2m ;fi
if [ $ivar == 6 ] ; then var=ShMF ;fi
if [ $ivar == 7 ] ; then var=CgMF ;fi

ireg=1
nreg=3
while [ $ireg -le $nreg ]
do
 if [ $ireg == 1 ] ; then reg=SAM ;fi
 if [ $ireg == 2 ] ; then reg=RSU ;fi
 if [ $ireg == 3 ] ; then reg=GLB ;fi

 echo $reg $var
 fig1=${reg}_${var}
 export fcst_files=`/bin/ls -1tr ${fig1}*.png`
 if [[ ${noutputs} -lt 7 ]] ; then 
    #${DIR_imagemagick}/montage -tile 3x2 -geometry +1+1 -background white $fcst_files  ${fig1}.x.png
     montage -tile 3x2 -geometry +1+1 -background white $fcst_files  ${fig1}.x.png

 fi
 if [[ ${noutputs} -gt 6 ]] ; then 
   #${DIR_imagemagick}/montage -tile 3x3 -geometry +1+1 -background white $fcst_files  ${fig1}.x.png
                       montage -tile 3x3 -geometry +1+1 -background white $fcst_files  ${fig1}.x.png
 fi

 ireg=`expr $ireg + 1`
done
ivar=`expr $ivar + 1`
done


#${DIR_imagemagick}/montage  -tile 3x2 -geometry +1+1 -background white \
#SAM_NCV_Prec.png SAM_CNV_Prec.png SAM_TOT_Prec.png \
#GLB_NCV_Prec.png GLB_CNV_Prec.png GLB_TOT_Prec.png \
#tar_precips.png 


/usr/bin/tar -cvf ${expname}_${start_data_dir}.tar *.x*.png mean_precips_*png
/usr/bin/mv ${expname}_${start_data_dir}.tar /home/saulo.freitas
/bin/mv *.png pngs



###############################################################
exit 
