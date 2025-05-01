#!/bin/bash
if [ -z "$1" ] ;then
    echo enter with expname for make_gf_nml
    exit
fi

exp=${1}
#  ---------------------------------------------------------------------------------------------------
if [ $exp == GFdef_1.4 ] ; then 
       /bin/cp GF_ConvPar_nml.bash GF_ConvPar_nml_${exp}.bash ; chmod +x  GF_ConvPar_nml_${exp}.bash
       
       exit
fi
#  ---------------------------------------------------------------------------------------------------
if [ $exp == GFdef_1.4_TcpL1h ] ; then 
       /bin/cp GF_ConvPar_nml.bash GF_ConvPar_nml_${exp}.bash ; chmod +x  GF_ConvPar_nml_${exp}.bash
       #sed -i  's/cum_HEI_UPDF_LAND = 0.55, 0.1, 0.55,/cum_HEI_UPDF_LAND = 0.65, 0.1, 0.55,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/cum_HEI_UPDF_OCEAN= 0.55, 0.1, 0.55,/cum_HEI_UPDF_OCEAN= 0.65, 0.1, 0.55,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/n_cldrop          = 100.,/n_cldrop          = 50.,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/overshoot         = 0.0,/overshoot         = 0.10,/g' GF_ConvPar_nml_${exp}.bash
       sed -i  's/tau_land_cp       = 1800.,/tau_land_cp       = 3600.,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/cum_entr_rate   = 0.8e-3, 3.e-3, 1.e-3,/cum_entr_rate   = 0.6e-3, 2.e-3, 1.e-3,/g' GF_ConvPar_nml_${exp}.bash    
       #sed -i  's/overshoot         = 0.0,/overshoot         = 0.10,/g' GF_ConvPar_nml_${exp}.bash

       exit
fi


#  ---------------------------------------------------------------------------------------------------
#if [ $exp == GFberyl_HDnO40L40_Smo2_OvSh20_sgsw2_TS30_ENT0.8M3_Tcp30mn_sigf0.65_edtL40 ] ; then 
if [ $exp == GFstart ] ; then 

/bin/cp GF_ConvPar_nml.bash GF_ConvPar_nml_${exp}.bash ; chmod +x  GF_ConvPar_nml_${exp}.bash
        #sed -i  's/cum_MAX_EDT_OCEAN = 0.5, 0.0, 0.5,/cum_MAX_EDT_OCEAN = 0.4, 0.0, 0.4,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/cum_MAX_EDT_LAND  = 0.5, 0.0, 0.2,/cum_MAX_EDT_LAND  = 0.4, 0.0, 0.8,/g' GF_ConvPar_nml_${exp}.bash
        #sed -i  's/dicycle           = 1,/dicycle           = 0,/g' GF_ConvPar_nml_${exp}.bash    
        sed -i  's/cum_t_star        = 4., -99., -99.,/cum_t_star = 30., -99., -99.,/g' GF_ConvPar_nml_${exp}.bash 
sed -i  's/cum_entr_rate     = 6.0e-4, 3.e-3, 1.e-3,/cum_entr_rate   = 0.8e-3, 2.e-3, 1.e-3,/g' GF_ConvPar_nml_${exp}.bash    
#sed -i  's/cum_entr_rate     = 6.0e-4, 3.e-3, 1.e-3,/cum_entr_rate   = 1.2e-3, 2.e-3, 1.e-3,/g' GF_ConvPar_nml_${exp}.bash    

        sed -i  's/icumulus_gf       = 1,1,2,/icumulus_gf   = 1,1,0,/g' GF_ConvPar_nml_${exp}.bash
        #sed -i  's/overshoot           = 0.20,/overshoot         = 0.0,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/cum_ave_layer     = 50.,30.,30.,/cum_ave_layer     = 30.,15.,30.,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/cum_min_cloud_depth = 200.,70.,150.,/cum_min_cloud_depth = 200.,50.,150.,/g' GF_ConvPar_nml_${exp}.bash    
        sed -i  's/cum_HEI_DOWN_LAND = 0.60, 0.0, 0.35,/cum_HEI_DOWN_LAND = 0.40, 0.0, 0.35,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/cum_HEI_DOWN_OCEAN= 0.50, 0.0, 0.35,/cum_HEI_DOWN_OCEAN= 0.40, 0.0, 0.35,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/max_tq_tend       = 300.,/max_tq_tend       = 200.,/g' GF_ConvPar_nml_${exp}.bash    
        #sed -i  's/tau_mid           = 3600.,/tau_mid = 3600.,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/cum_HEI_UPDF_LAND = 0.55, 0.1, 0.55,/cum_HEI_UPDF_LAND = 0.40, 0.1, 0.55,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/cum_HEI_UPDF_OCEAN= 0.55, 0.1, 0.55,/cum_HEI_UPDF_OCEAN= 0.40, 0.1, 0.55,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/n_cldrop          = 100.,/n_cldrop          = 50.,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/autoconv          = 4,/autoconv          = 1,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/qrc_crit          = 6.0e-4,/qrc_crit          = 1.0e-3,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/sgs_w_timescale   = 1,/sgs_w_timescale   = 2,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/cum_use_excess    = 1,1,1,/cum_use_excess    = 2,1,1,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/use_rhu_ctrl_entr = 1,/use_rhu_ctrl_entr = 0,/g' GF_ConvPar_nml_${exp}.bash    
       #sed -i  's/cum_use_smooth_tend = 2,2,2,/cum_use_smooth_tend   = 1,1,1,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/n_cldrop          = 100.,/n_cldrop          = 150.,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/use_linear_subcl_mf = 1,/use_linear_subcl_mf = 0,/g' GF_ConvPar_nml_${exp}.bash  
       #sed -i  's/c0_shal           = 0.0e-3,/c0_shal           = 1.0e-3,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/tau_ocea_cp       = 7200.,/tau_ocea_cp       = 3600.,/g' GF_ConvPar_nml_${exp}.bash
       sed -i  's/tau_land_cp       = 7200.,/tau_land_cp       = 1800.,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/use_scale_dep     = 1,/use_scale_dep     = 0,/g' GF_ConvPar_nml_${exp}.bash 
        sed -i  's/sig_factor        = 0.22,/sig_factor        = 0.65,/g' GF_ConvPar_nml_${exp}.bash         

exit
fi

#  ---------------------------------------------------------------------------------------------------
#if [ $exp == GFstart_shAL30CD70E3M3c01M3_OvShOFF ] ; then 
#if [ $exp == GFstart2_qrcL1M4O3M4c01M3_TcpO1h ] ; then 
if [ $exp == GFstart3_EdtO60L30 ] ; then 

/bin/cp GF_ConvPar_nml.bash GF_ConvPar_nml_${exp}.bash ; chmod +x  GF_ConvPar_nml_${exp}.bash
        sed -i  's/cum_MAX_EDT_OCEAN = 0.5, 0.0, 0.5,/cum_MAX_EDT_OCEAN = 0.6, 0.0, 0.4,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/cum_MAX_EDT_LAND  = 0.5, 0.0, 0.2,/cum_MAX_EDT_LAND  = 0.3, 0.0, 0.8,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/dicycle           = 1,/dicycle           = 0,/g' GF_ConvPar_nml_${exp}.bash    
        sed -i  's/cum_t_star        = 4., -99., -99.,/cum_t_star = 30., -99., -99.,/g' GF_ConvPar_nml_${exp}.bash 
        sed -i  's/cum_entr_rate     = 6.0e-4, 3.e-3, 1.e-3,/cum_entr_rate   = 0.8e-3, 3.e-3, 1.e-3,/g' GF_ConvPar_nml_${exp}.bash    

        sed -i  's/icumulus_gf       = 1,1,2,/icumulus_gf   = 1,1,0,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/overshoot           = 0.20,/overshoot         = 0.0,/g' GF_ConvPar_nml_${exp}.bash

#       sed -i  's/cum_ave_layer     = 50.,30.,30.,/cum_ave_layer     = 30.,15.,30.,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/cum_ave_layer     = 50.,30.,30.,/cum_ave_layer     = 30.,30.,30.,/g' GF_ConvPar_nml_${exp}.bash
#       sed -i  's/cum_min_cloud_depth = 200.,70.,150.,/cum_min_cloud_depth = 200.,50.,150.,/g' GF_ConvPar_nml_${exp}.bash    
        sed -i  's/cum_min_cloud_depth = 200.,70.,150.,/cum_min_cloud_depth = 200.,70.,150.,/g' GF_ConvPar_nml_${exp}.bash    

        sed -i  's/cum_HEI_DOWN_LAND = 0.60, 0.0, 0.35,/cum_HEI_DOWN_LAND = 0.40, 0.0, 0.35,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/cum_HEI_DOWN_OCEAN= 0.50, 0.0, 0.35,/cum_HEI_DOWN_OCEAN= 0.40, 0.0, 0.35,/g' GF_ConvPar_nml_${exp}.bash

        sed -i  's/max_tq_tend       = 300.,/max_tq_tend       = 200.,/g' GF_ConvPar_nml_${exp}.bash    

       #sed -i  's/tau_mid           = 3600.,/tau_mid = 3600.,/g' GF_ConvPar_nml_${exp}.bash

       #sed -i  's/cum_HEI_UPDF_LAND = 0.55, 0.1, 0.55,/cum_HEI_UPDF_LAND = 0.40, 0.1, 0.55,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/cum_HEI_UPDF_OCEAN= 0.55, 0.1, 0.55,/cum_HEI_UPDF_OCEAN= 0.40, 0.1, 0.55,/g' GF_ConvPar_nml_${exp}.bash

       #sed -i  's/n_cldrop          = 100.,/n_cldrop          = 50.,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/autoconv          = 4,/autoconv          = 1,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/qrc_crit          = 1.0e-3,/qrc_crit          = 1.0e-4,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/sgs_w_timescale   = 1,/sgs_w_timescale   = 2,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/cum_use_excess    = 1,1,1,/cum_use_excess    = 2,1,1,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/use_rhu_ctrl_entr = 1,/use_rhu_ctrl_entr = 0,/g' GF_ConvPar_nml_${exp}.bash    
       #sed -i  's/cum_use_smooth_tend = 2,2,2,/cum_use_smooth_tend   = 1,1,1,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/n_cldrop          = 100.,/n_cldrop          = 150.,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/use_linear_subcl_mf = 1,/use_linear_subcl_mf = 0,/g' GF_ConvPar_nml_${exp}.bash  
        sed -i  's/c0_shal           = 0.0e-3,/c0_shal           = 1.0e-3,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/tau_ocea_cp       = 7200.,/tau_ocea_cp       = 3600.,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/tau_land_cp       = 7200.,/tau_land_cp       = 1800.,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/use_scale_dep     = 1,/use_scale_dep     = 0,/g' GF_ConvPar_nml_${exp}.bash 
        sed -i  's/sig_factor        = 0.22,/sig_factor        = 0.65,/g' GF_ConvPar_nml_${exp}.bash         

exit
fi

#  ---------------------------------------------------------------------------------------------------
#if [ $exp == GFstart_shAL30CD70E3M3c01M3_OvShOFF ] ; then 
if [ $exp == GFstart2 ] ; then 

/bin/cp GF_ConvPar_nml.bash GF_ConvPar_nml_${exp}.bash ; chmod +x  GF_ConvPar_nml_${exp}.bash
       #sed -i  's/cum_MAX_EDT_OCEAN = 0.5, 0.0, 0.5,/cum_MAX_EDT_OCEAN = 0.4, 0.0, 0.4,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/cum_MAX_EDT_LAND  = 0.5, 0.0, 0.2,/cum_MAX_EDT_LAND  = 0.4, 0.0, 0.8,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/dicycle           = 1,/dicycle           = 0,/g' GF_ConvPar_nml_${exp}.bash    
        sed -i  's/cum_t_star        = 4., -99., -99.,/cum_t_star = 30., -99., -99.,/g' GF_ConvPar_nml_${exp}.bash 
        sed -i  's/cum_entr_rate     = 6.0e-4, 3.e-3, 1.e-3,/cum_entr_rate   = 0.8e-3, 3.e-3, 1.e-3,/g' GF_ConvPar_nml_${exp}.bash    

        sed -i  's/icumulus_gf       = 1,1,2,/icumulus_gf   = 1,1,0,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/overshoot           = 0.20,/overshoot         = 0.0,/g' GF_ConvPar_nml_${exp}.bash

#       sed -i  's/cum_ave_layer     = 50.,30.,30.,/cum_ave_layer     = 30.,15.,30.,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/cum_ave_layer     = 50.,30.,30.,/cum_ave_layer     = 30.,30.,30.,/g' GF_ConvPar_nml_${exp}.bash
#       sed -i  's/cum_min_cloud_depth = 200.,70.,150.,/cum_min_cloud_depth = 200.,50.,150.,/g' GF_ConvPar_nml_${exp}.bash    
        sed -i  's/cum_min_cloud_depth = 200.,70.,150.,/cum_min_cloud_depth = 200.,70.,150.,/g' GF_ConvPar_nml_${exp}.bash    

        sed -i  's/cum_HEI_DOWN_LAND = 0.60, 0.0, 0.35,/cum_HEI_DOWN_LAND = 0.40, 0.0, 0.35,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/cum_HEI_DOWN_OCEAN= 0.50, 0.0, 0.35,/cum_HEI_DOWN_OCEAN= 0.40, 0.0, 0.35,/g' GF_ConvPar_nml_${exp}.bash

        sed -i  's/max_tq_tend       = 300.,/max_tq_tend       = 200.,/g' GF_ConvPar_nml_${exp}.bash    

       #sed -i  's/tau_mid           = 3600.,/tau_mid = 3600.,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/cum_HEI_UPDF_LAND = 0.55, 0.1, 0.55,/cum_HEI_UPDF_LAND = 0.40, 0.1, 0.55,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/cum_HEI_UPDF_OCEAN= 0.55, 0.1, 0.55,/cum_HEI_UPDF_OCEAN= 0.40, 0.1, 0.55,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/n_cldrop          = 100.,/n_cldrop          = 50.,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/autoconv          = 4,/autoconv          = 1,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/qrc_crit          = 6.0e-4,/qrc_crit          = 1.0e-3,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/sgs_w_timescale   = 1,/sgs_w_timescale   = 2,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/cum_use_excess    = 1,1,1,/cum_use_excess    = 2,1,1,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/use_rhu_ctrl_entr = 1,/use_rhu_ctrl_entr = 0,/g' GF_ConvPar_nml_${exp}.bash    
       #sed -i  's/cum_use_smooth_tend = 2,2,2,/cum_use_smooth_tend   = 1,1,1,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/n_cldrop          = 100.,/n_cldrop          = 150.,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/use_linear_subcl_mf = 1,/use_linear_subcl_mf = 0,/g' GF_ConvPar_nml_${exp}.bash  
        sed -i  's/c0_shal           = 0.0e-3,/c0_shal           = 1.0e-3,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/tau_ocea_cp       = 7200.,/tau_ocea_cp       = 3600.,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/tau_land_cp       = 7200.,/tau_land_cp       = 1800.,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/use_scale_dep     = 1,/use_scale_dep     = 0,/g' GF_ConvPar_nml_${exp}.bash 
        sed -i  's/sig_factor        = 0.22,/sig_factor        = 0.65,/g' GF_ConvPar_nml_${exp}.bash         

exit
fi



#  ---------------------------------------------------------------------------------------------------
#if [ $exp == GFberyl_HDnO40L40_Smo2_OvSh20_sgsw2_TS30_ENT0.8M3_Tcp30mn_sigf0.65_edtL40 ] ; then 
if [ $exp == GFstart_EdtO60HdnO35 ] ; then 

/bin/cp GF_ConvPar_nml.bash GF_ConvPar_nml_${exp}.bash ; chmod +x  GF_ConvPar_nml_${exp}.bash
        #sed -i  's/cum_MAX_EDT_OCEAN = 0.5, 0.0, 0.5,/cum_MAX_EDT_OCEAN = 0.4, 0.0, 0.4,/g' GF_ConvPar_nml_${exp}.bash
         sed -i  's/cum_MAX_EDT_OCEAN = 0.5, 0.0, 0.5,/cum_MAX_EDT_OCEAN = 0.6, 0.0, 0.4,/g' GF_ConvPar_nml_${exp}.bash

        sed -i  's/cum_MAX_EDT_LAND  = 0.5, 0.0, 0.2,/cum_MAX_EDT_LAND  = 0.4, 0.0, 0.8,/g' GF_ConvPar_nml_${exp}.bash
        #sed -i  's/dicycle           = 1,/dicycle           = 0,/g' GF_ConvPar_nml_${exp}.bash    
        sed -i  's/cum_t_star        = 4., -99., -99.,/cum_t_star = 30., -99., -99.,/g' GF_ConvPar_nml_${exp}.bash 
sed -i  's/cum_entr_rate     = 6.0e-4, 3.e-3, 1.e-3,/cum_entr_rate   = 0.8e-3, 2.e-3, 1.e-3,/g' GF_ConvPar_nml_${exp}.bash    
#sed -i  's/cum_entr_rate     = 6.0e-4, 3.e-3, 1.e-3,/cum_entr_rate   = 1.2e-3, 2.e-3, 1.e-3,/g' GF_ConvPar_nml_${exp}.bash    
        sed -i  's/icumulus_gf       = 1,1,2,/icumulus_gf   = 1,1,0,/g' GF_ConvPar_nml_${exp}.bash
        #sed -i  's/overshoot           = 0.20,/overshoot         = 0.0,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/cum_ave_layer     = 50.,30.,30.,/cum_ave_layer     = 30.,15.,30.,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/cum_min_cloud_depth = 200.,70.,150.,/cum_min_cloud_depth = 200.,50.,150.,/g' GF_ConvPar_nml_${exp}.bash    
        sed -i  's/cum_HEI_DOWN_LAND = 0.60, 0.0, 0.35,/cum_HEI_DOWN_LAND = 0.40, 0.0, 0.35,/g' GF_ConvPar_nml_${exp}.bash

#srf    sed -i  's/cum_HEI_DOWN_OCEAN= 0.50, 0.0, 0.35,/cum_HEI_DOWN_OCEAN= 0.40, 0.0, 0.35,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/cum_HEI_DOWN_OCEAN= 0.50, 0.0, 0.35,/cum_HEI_DOWN_OCEAN= 0.35, 0.0, 0.35,/g' GF_ConvPar_nml_${exp}.bash

        sed -i  's/max_tq_tend       = 300.,/max_tq_tend       = 200.,/g' GF_ConvPar_nml_${exp}.bash    
        #sed -i  's/tau_mid           = 3600.,/tau_mid = 3600.,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/cum_HEI_UPDF_LAND = 0.55, 0.1, 0.55,/cum_HEI_UPDF_LAND = 0.40, 0.1, 0.55,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/cum_HEI_UPDF_OCEAN= 0.55, 0.1, 0.55,/cum_HEI_UPDF_OCEAN= 0.40, 0.1, 0.55,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/n_cldrop          = 100.,/n_cldrop          = 50.,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/autoconv          = 4,/autoconv          = 1,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/qrc_crit          = 6.0e-4,/qrc_crit          = 1.0e-3,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/sgs_w_timescale   = 1,/sgs_w_timescale   = 2,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/cum_use_excess    = 1,1,1,/cum_use_excess    = 2,1,1,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/use_rhu_ctrl_entr = 1,/use_rhu_ctrl_entr = 0,/g' GF_ConvPar_nml_${exp}.bash    
       #sed -i  's/cum_use_smooth_tend = 2,2,2,/cum_use_smooth_tend   = 1,1,1,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/n_cldrop          = 100.,/n_cldrop          = 150.,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/use_linear_subcl_mf = 1,/use_linear_subcl_mf = 0,/g' GF_ConvPar_nml_${exp}.bash  
       #sed -i  's/c0_shal           = 0.0e-3,/c0_shal           = 1.0e-3,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/tau_ocea_cp       = 7200.,/tau_ocea_cp       = 3600.,/g' GF_ConvPar_nml_${exp}.bash
       sed -i  's/tau_land_cp       = 7200.,/tau_land_cp       = 1800.,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/use_scale_dep     = 1,/use_scale_dep     = 0,/g' GF_ConvPar_nml_${exp}.bash 
        sed -i  's/sig_factor        = 0.22,/sig_factor        = 0.65,/g' GF_ConvPar_nml_${exp}.bash         

exit
fi







#  ---------------------------------------------------------------------------------------------------
if [ $exp == GFberyl_HDnO40L40_Smo2_OvSh20_sgsw2_TS10 ] ; then 

/bin/cp GF_ConvPar_nml.bash GF_ConvPar_nml_${exp}.bash ; chmod +x  GF_ConvPar_nml_${exp}.bash
        #sed -i  's/cum_MAX_EDT_OCEAN = 0.5, 0.0, 0.5,/cum_MAX_EDT_OCEAN = 0.4, 0.0, 0.4,/g' GF_ConvPar_nml_${exp}.bash
        #sed -i  's/cum_MAX_EDT_LAND  = 0.5, 0.0, 0.2,/cum_MAX_EDT_LAND  = 0.6, 0.0, 0.8,/g' GF_ConvPar_nml_${exp}.bash
        #sed -i  's/dicycle           = 1,/dicycle           = 0,/g' GF_ConvPar_nml_${exp}.bash    
        sed -i  's/cum_t_star        = 4., -99., -99.,/cum_t_star = 10., -99., -99.,/g' GF_ConvPar_nml_${exp}.bash 
        sed -i  's/cum_entr_rate     = 6.0e-4, 3.e-3, 1.e-3,/cum_entr_rate   = 6.0e-4, 2.e-3, 1.e-3,/g' GF_ConvPar_nml_${exp}.bash    
       
        sed -i  's/icumulus_gf       = 1,1,2,/icumulus_gf   = 1,1,0,/g' GF_ConvPar_nml_${exp}.bash
        #sed -i  's/overshoot           = 0.20,/overshoot         = 0.0,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/cum_ave_layer     = 50.,30.,30.,/cum_ave_layer     = 30.,15.,30.,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/cum_min_cloud_depth = 200.,70.,150.,/cum_min_cloud_depth = 200.,50.,150.,/g' GF_ConvPar_nml_${exp}.bash    
        sed -i  's/cum_HEI_DOWN_LAND = 0.60, 0.0, 0.35,/cum_HEI_DOWN_LAND = 0.40, 0.0, 0.35,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/cum_HEI_DOWN_OCEAN= 0.50, 0.0, 0.35,/cum_HEI_DOWN_OCEAN= 0.40, 0.0, 0.35,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/max_tq_tend       = 300.,/max_tq_tend       = 200.,/g' GF_ConvPar_nml_${exp}.bash    
        #sed -i  's/tau_mid           = 3600.,/tau_mid = 3600.,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/cum_HEI_UPDF_LAND = 0.55, 0.1, 0.55,/cum_HEI_UPDF_LAND = 0.40, 0.1, 0.55,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/cum_HEI_UPDF_OCEAN= 0.55, 0.1, 0.55,/cum_HEI_UPDF_OCEAN= 0.40, 0.1, 0.55,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/n_cldrop          = 100.,/n_cldrop          = 50.,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/autoconv          = 4,/autoconv          = 1,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/qrc_crit          = 6.0e-4,/qrc_crit          = 1.0e-3,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/sgs_w_timescale   = 1,/sgs_w_timescale   = 2,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/cum_use_excess    = 1,1,1,/cum_use_excess    = 2,1,1,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/use_rhu_ctrl_entr = 1,/use_rhu_ctrl_entr = 0,/g' GF_ConvPar_nml_${exp}.bash    
       #sed -i  's/cum_use_smooth_tend = 2,2,2,/cum_use_smooth_tend   = 1,1,1,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/n_cldrop          = 100.,/n_cldrop          = 150.,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/use_linear_subcl_mf = 1,/use_linear_subcl_mf = 0,/g' GF_ConvPar_nml_${exp}.bash  
       #sed -i  's/c0_shal           = 0.0e-3,/c0_shal           = 1.0e-3,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/tau_ocea_cp       = 7200.,/tau_ocea_cp       = 3600.,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/tau_land_cp       = 7200.,/tau_land_cp       = 2400.,/g' GF_ConvPar_nml_${exp}.bash
 
exit
fi

#  ---------------------------------------------------------------------------------------------------
if [ $exp == GFdef_target_Tcp3h ] ; then 

/bin/cp GF_ConvPar_nml.bash GF_ConvPar_nml_${exp}.bash ; chmod +x  GF_ConvPar_nml_${exp}.bash
        #sed -i  's/cum_MAX_EDT_OCEAN = 0.5, 0.0, 0.5,/cum_MAX_EDT_OCEAN = 0.4, 0.0, 0.4,/g' GF_ConvPar_nml_${exp}.bash
        #sed -i  's/cum_MAX_EDT_LAND  = 0.5, 0.0, 0.2,/cum_MAX_EDT_LAND  = 0.6, 0.0, 0.8,/g' GF_ConvPar_nml_${exp}.bash

        #sed -i  's/cum_entr_rate     = 6.0e-4, 3.e-3, 1.e-3,/cum_entr_rate   = 6.0e-4, 2.e-3, 1.e-3,/g' GF_ConvPar_nml_${exp}.bash    
        sed -i  's/icumulus_gf       = 1,1,2,/icumulus_gf   = 1,1,0,/g' GF_ConvPar_nml_${exp}.bash
        #sed -i  's/overshoot           = 0.20,/overshoot         = 0.0,/g' GF_ConvPar_nml_${exp}.bash
        #sed -i  's/cum_ave_layer     = 50.,30.,30.,/cum_ave_layer     = 30.,15.,30.,/g' GF_ConvPar_nml_${exp}.bash
        #sed -i  's/cum_min_cloud_depth = 200.,70.,150.,/cum_min_cloud_depth = 200.,50.,150.,/g' GF_ConvPar_nml_${exp}.bash    
        sed -i  's/cum_HEI_DOWN_LAND = 0.60, 0.0, 0.35,/cum_HEI_DOWN_LAND = 0.45, 0.0, 0.35,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/cum_HEI_DOWN_OCEAN= 0.50, 0.0, 0.35,/cum_HEI_DOWN_OCEAN= 0.45, 0.0, 0.35,/g' GF_ConvPar_nml_${exp}.bash
        #sed -i  's/max_tq_tend       = 300.,/max_tq_tend       = 200.,/g' GF_ConvPar_nml_${exp}.bash    
       sed -i  's/tau_ocea_cp       = 7200.,/tau_ocea_cp       = 10800.,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/tau_land_cp       = 7200.,/tau_land_cp       = 3600.,/g' GF_ConvPar_nml_${exp}.bash
 
exit
fi
#  ---
#  ---------------------------------------------------------------------------------------------------
if [ $exp == GFberyl_avely20 ] ; then 

/bin/cp GF_ConvPar_nml.bash GF_ConvPar_nml_${exp}.bash ; chmod +x  GF_ConvPar_nml_${exp}.bash
        sed -i  's/cum_MAX_EDT_OCEAN = 0.5, 0.0, 0.5,/cum_MAX_EDT_OCEAN = 0.5, 0.0, 0.4,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/cum_MAX_EDT_LAND  = 0.5, 0.0, 0.2,/cum_MAX_EDT_LAND  = 0.5, 0.0, 0.8,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/cum_entr_rate     = 6.0e-4, 3.e-3, 1.e-3,/cum_entr_rate   = 6.0e-4, 2.e-3, 1.e-3,/g' GF_ConvPar_nml_${exp}.bash    
        sed -i  's/icumulus_gf       = 1,1,2,/icumulus_gf   = 1,1,0,/g' GF_ConvPar_nml_${exp}.bash     
        sed -i  's/cum_use_smooth_tend = 2,2,2,/cum_use_smooth_tend   = 1,1,1,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/overshoot           = 0.20,/overshoot         = 0.0,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/cum_ave_layer     = 50.,30.,30.,/cum_ave_layer     = 20.,15.,30.,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/cum_min_cloud_depth = 200.,70.,150.,/cum_min_cloud_depth = 200.,50.,150.,/g' GF_ConvPar_nml_${exp}.bash    
        sed -i  's/cum_HEI_DOWN_LAND = 0.60, 0.0, 0.35,/cum_HEI_DOWN_LAND = 0.40, 0.0, 0.35,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/cum_HEI_DOWN_OCEAN= 0.50, 0.0, 0.35,/cum_HEI_DOWN_OCEAN= 0.35, 0.0, 0.35,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/max_tq_tend       = 300.,/max_tq_tend       = 200.,/g' GF_ConvPar_nml_${exp}.bash    
        sed -i  's/tau_mid           = 3600.,/tau_mid = 1200.,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/cum_HEI_UPDF_LAND = 0.55, 0.1, 0.55,/cum_HEI_UPDF_LAND = 0.55, 0.1, 0.55,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/cum_HEI_UPDF_OCEAN= 0.55, 0.1, 0.55,/cum_HEI_UPDF_OCEAN= 0.55, 0.1, 0.55,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/n_cldrop          = 100.,/n_cldrop          = 50.,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/autoconv          = 4,/autoconv          = 1,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/qrc_crit          = 6.0e-4,/qrc_crit          = 1.0e-3,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/sgs_w_timescale   = 0,/sgs_w_timescale   = 1,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/cum_use_excess    = 1,1,1,/cum_use_excess    = 2,1,1,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/use_rhu_ctrl_entr = 1,/use_rhu_ctrl_entr = 0,/g' GF_ConvPar_nml_${exp}.bash    
       #sed -i  's/n_cldrop          = 100.,/n_cldrop          = 150.,/g' GF_ConvPar_nml_${exp}.bash
       # sed -i  's/tau_ocea_cp       = 7200.,/tau_ocea_cp       = 10800.,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/tau_land_cp       = 7200.,/tau_land_cp       = 3600.,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/use_linear_subcl_mf = 1,/use_linear_subcl_mf = 0,/g' GF_ConvPar_nml_${exp}.bash  
       #sed -i  's/c0_shal           = 0.0e-3,/c0_shal           = 1.0e-3,/g' GF_ConvPar_nml_${exp}.bash
exit
fi



#  ---------------------------------------------------------------------------------------------------
if [ $exp == GFprp_SH1_En3Shlay30dp70_DpHup55_Cg2T1hEDT50Hdn35REF_TQ300  ] ; then 

/bin/cp GF_ConvPar_nml.bash GF_ConvPar_nml_${exp}.bash ; chmod +x  GF_ConvPar_nml_${exp}.bash
        sed -i  's/cum_entr_rate     = 6.0e-4, 3.e-3, 1.e-3,/cum_entr_rate   = 6.0e-4, 3.e-3, 1.e-3,/g' GF_ConvPar_nml_${exp}.bash    
        sed -i  's/icumulus_gf       = 1,1,2,/icumulus_gf   = 1,1,2,/g' GF_ConvPar_nml_${exp}.bash     

        sed -i  's/overshoot           = 0.20,/overshoot         = 0.2,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/cum_ave_layer     = 50.,30.,30.,/cum_ave_layer     = 50.,30.,30.,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/cum_min_cloud_depth = 200.,70.,150.,/cum_min_cloud_depth = 200.,70.,150.,/g' GF_ConvPar_nml_${exp}.bash    
        sed -i  's/cum_HEI_DOWN_LAND = 0.60, 0.0, 0.35,/cum_HEI_DOWN_LAND = 0.60, 0.0, 0.35,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/cum_HEI_DOWN_OCEAN= 0.50, 0.0, 0.35,/cum_HEI_DOWN_OCEAN= 0.50, 0.0, 0.35,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/max_tq_tend       = 300.,/max_tq_tend       = 300.,/g' GF_ConvPar_nml_${exp}.bash    
        sed -i  's/tau_mid           = 3600.,/tau_mid = 3600.,/g' GF_ConvPar_nml_${exp}.bash
        sed -i  's/cum_MAX_EDT_OCEAN = 0.5, 0.0, 0.5,/cum_MAX_EDT_OCEAN = 0.5, 0.0, 0.5,/g' GF_ConvPar_nml_${exp}.bash

       #sed -i  's/cum_MAX_EDT_LAND  = 0.5, 0.0, 0.2,/cum_MAX_EDT_LAND  = 0.5, 0.0, 0.5,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/cum_HEI_UPDF_LAND = 0.55, 0.1, 0.55,/cum_HEI_UPDF_LAND = 0.40, 0.1, 0.55,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/cum_HEI_UPDF_OCEAN= 0.55, 0.1, 0.55,/cum_HEI_UPDF_OCEAN= 0.40, 0.1, 0.55,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/n_cldrop          = 100.,/n_cldrop          = 50.,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/autoconv          = 4,/autoconv          = 1,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/qrc_crit          = 6.0e-4,/qrc_crit          = 1.0e-3,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/sgs_w_timescale   = 0,/sgs_w_timescale   = 1,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/cum_use_excess    = 1,1,1,/cum_use_excess    = 2,1,1,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/use_rhu_ctrl_entr = 1,/use_rhu_ctrl_entr = 0,/g' GF_ConvPar_nml_${exp}.bash    
       #sed -i  's/cum_use_smooth_tend = 2,2,2,/cum_use_smooth_tend   = 2,2,2,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/n_cldrop          = 100.,/n_cldrop          = 150.,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/tau_ocea_cp       = 7200.,/tau_ocea_cp       = 3600.,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/tau_land_cp       = 7200.,/tau_land_cp       = 3600.,/g' GF_ConvPar_nml_${exp}.bash
       #sed -i  's/use_linear_subcl_mf = 1,/use_linear_subcl_mf = 0,/g' GF_ConvPar_nml_${exp}.bash  
       #sed -i  's/c0_shal           = 0.0e-3,/c0_shal           = 1.0e-3,/g' GF_ConvPar_nml_${exp}.bash
exit
fi

#  ---------------------------------------------------------------------------------------------------

#
# backup of options - do not delete
# ----------------------------------------------------------------------------------------------------
exit 



#         sed -i 's/convection_tracer = 0,/convection_tracer = 1,/g' GF_ConvPar_nml_${exp}.bash    
#         sed -i  's/use_rhu_ctrl_entr = 0,/use_rhu_ctrl_entr = 0,/g' GF_ConvPar_nml_${exp}.bash    
#         sed -i  's/cum_entr_rate     = 6.0e-4, 2.e-3, 1.e-3,/cum_entr_rate   = 1.0e-4, 1.e-3, 2.e-3, /g' GF_ConvPar_nml_${exp}.bash    
#         sed -i  's/cum_min_cloud_depth = 200.,100.,150.,/cum_min_cloud_depth = 0.,0.,0.,/g' GF_ConvPar_nml_${exp}.bash    
#         sed -i  's/add_coldpool_trig = 0,/add_coldpool_trig = 2,/g' GF_ConvPar_nml_${exp}.bash     
#         sed -i  's/add_coldpool_clos = 0,/add_coldpool_clos = 2,/g' GF_ConvPar_nml_${exp}.bash     
#         sed -i  's/use_pass_cloudvol = 0,/use_pass_cloudvol = 2,/g' GF_ConvPar_nml_${exp}.bash    
#         sed -i  's/use_memory        = 0,/use_memory        = 22,/g' GF_ConvPar_nml_${exp}.bash    
#         sed -i  's/use_lcl_ctrl_entr = 0,/use_lcl_ctrl_entr = 2,/g' GF_ConvPar_nml_${exp}.bash    
#         sed -i  's/cum_MAX_EDT_LAND  = 0.5, 0.0, 0.2,/cum_MAX_EDT_LAND  = 0.3, 0.0, 0.2,/g' GF_ConvPar_nml_${exp}.bash
#         sed -i  's/cum_MAX_EDT_OCEAN = 0.5, 0.0, 0.2,/cum_MAX_EDT_OCEAN = 0.3, 0.0, 0.2,/g' GF_ConvPar_nml_${exp}.bash
#         sed -i  's/dicycle           = 0,/dicycle           = 1,/g' GF_ConvPar_nml_${exp}.bash    

        # sed -i  's/c0_deep           = 1.0e-3,/c0_deep      = 2.0e-3,/g' GF_ConvPar_nml_${exp}.bash    
        # sed -i  's/output_sound      = 0,/output_sound      = 3,/g' GF_ConvPar_nml_${exp}.bash    
        # sed -i  's/closure_choice    = 10,10,3,/closure_choice = 10,1,3,/g' GF_ConvPar_nml_${exp}.bash    
        # sed -i  's/icumulus_gf       = 1,1,0,/icumulus_gf   = 1,1,1,/g' GF_ConvPar_nml_${exp}.bash     
        # sed -i  's/downdraft         = 1,/downdraft         = 0,/g' GF_ConvPar_nml_${exp}.bash    
        # sed -i  's/use_rebcb         = 1,/use_rebcb         = 0,/g' GF_ConvPar_nml_${exp}.bash    
#        sed -i  's/convection_tracer = 1,/convection_tracer = 0,/g' GF_ConvPar_nml_${exp}.bash    
#        sed -i  's/add_coldpool_trig = 2,/add_coldpool_trig = 0,/g' GF_ConvPar_nml_${exp}.bash     
#        sed -i  's/add_coldpool_clos = 2,/add_coldpool_clos = 02,/g' GF_ConvPar_nml_${exp}.bash     
#        sed -i  's/use_memory        = 222,/use_memory      = 0,/g' GF_ConvPar_nml_${exp}.bash    

