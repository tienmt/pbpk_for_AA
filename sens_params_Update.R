# 2. Define a character vector of ONLY the parameters that actually exist in the model
sens_params <- c(
  'pAA_TB', 'pAA_LiB', 'pAA_KiB'  ,
  'k_onAA_T', 'k_onAA_B', 'k_onAA_Li'  ,
  'pAA_TB', 'pAA_LiB' ,  'pAA_KiB'   ,
  'k_onAA_T', 'k_onAA_B' , 'k_onAA_Li', 'k_onAA_Ki'   , 
  'k_onAA_Ki',  'k_exc_AAMA', 'V_max_EH', 'KM_EH'
  ,
  'K_FORM_AA_VAL', 'K_REM_AA_VAL'  ,  'V_max_EH' 
  , 
  'k_cl_GSH', 'V_max_p450', 'KM_p450'
  ,
  'pGA_TB', 'pGA_LiB', 'pGA_KiB'  ,
  'k_onGA_T', 'k_onGA_B', 'k_onGA_Li'  ,
  'pGA_TB', 'pGA_LiB' ,  'pGA_KiB'    ,
  'k_onGA_T', 'k_onGA_B' , 'k_onGA_Li', 'k_onGA_Ki' 
)

##Here there are parameters that are mentioned two times and sometimes can lead 
## to incorrect sensitivity analysis. Also, I think that are missing few important
## parameters e.g. Vmax_AA_GSH for the AAMA sensitivity output, I suggest we change
#sens_params <- c(
#  'pAA_TB', 'pAA_LiB', 'pAA_KiB'  ,
#  'k_onAA_B', 'k_onAA_T', 'k_onAA_Li', 
#  'k_onAA_Ki',  'k_exc_AAMA', 'V_max_EH', 'KM_EH',
#  'K_FORM_AA_VAL', 'K_REM_AA_VAL'  ,
#  'k_cl_GSH', 'V_max_p450', 'KM_p450',
#  'pGA_TB', 'pGA_LiB', 'pGA_KiB'  ,
#  'k_onGA_T', 'k_onGA_B', 'k_onGA_Li', 'k_onGA_Ki',
#  'Vmax_AA_GSH', 'Vmax_GA_GSH','k_exc_GAMA', 'k_exc_GAOH',
#  'K_FORM_GA_VAL',  'K_REM_GA_VAL')   

