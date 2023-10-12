% This function returns the standard deviation of the optimal fibers
% The data come from Rajagopal et al. (2016). (Ward et al.)
% We used 0 when no data were available.
%
% Author: Israel Luis
% Date: 04/01/2023
%     
function data = getOptimalFiberSD(muscleNames)

    oflSD_data.glut_med1_r = 1.6;      oflSD_data.glmed1_r = 1.6;
    oflSD_data.glut_med2_r = 1.6;      oflSD_data.glmed2_r = 1.6;
    oflSD_data.glut_med3_r = 1.6;      oflSD_data.glmed3_r = 1.6;
    oflSD_data.glut_min1_r = NaN;      oflSD_data.glmin1_r = NaN;
    oflSD_data.glut_min2_r = NaN;      oflSD_data.glmin2_r = NaN;
    oflSD_data.glut_min3_r = NaN;      oflSD_data.glmin3_r = NaN;
    oflSD_data.semimem_r   = 1.8; 
    oflSD_data.semiten_r   = 4.1;   
    
    oflSD_data.bifemlh_r   = 2.6;      oflSD_data.bflh_r = 2.6;         oflSD_data.bflh140_r = 2.6; 
    oflSD_data.bifemsh_r   = 2.1;      oflSD_data.bfsh_r = 2.1;          oflSD_data.bfsh140_r = 2.1;  
    
    oflSD_data.sar_r       = 4.6;      oflSD_data.sart_r = 4.6;
    oflSD_data.add_mag1_r  = 3.4;      oflSD_data.addmagDist_r = 3.4;
    oflSD_data.add_mag2_r  = 3.0;      oflSD_data.addmagIsch_r = 3.0;
    oflSD_data.add_mag3_r  = 2.6;      oflSD_data.addmagMid_r  = 2.6;     
                                       oflSD_data.addmagProx_r = 2.0;
    oflSD_data.tfl_r       = NaN;
    oflSD_data.pect_r      = NaN;
    oflSD_data.grac_r      = 4.4;
    oflSD_data.glut_max1_r = 2.4;      oflSD_data.glmax1_r = 2.4;
    oflSD_data.glut_max2_r = 2.6;      oflSD_data.glmax2_r = 2.6;
    oflSD_data.glut_max3_r = 2.7;      oflSD_data.glmax3_r = 2.7;
    
    oflSD_data.iliacus_r   = 1.9;
    oflSD_data.psoas_r     = 1.7;
    oflSD_data.quad_fem_r  = NaN;
    oflSD_data.gem_r       = NaN;
    oflSD_data.peri_r      = NaN;      oflSD_data.piri_r = NaN;
    
    oflSD_data.rect_fem_r  = 1.3;      oflSD_data.recfem_r = 1.3;
    oflSD_data.vas_med_r = 2.3;        oflSD_data.vasmed_r = 2.3;
    oflSD_data.vas_int_r = 2.0;        oflSD_data.vasint_r = 2.0;
    oflSD_data.vas_lat_r = 1.8;        oflSD_data.vaslat_r = 1.8;        oflSD_data.vaslat140_r = 1.8; 
    
    oflSD_data.med_gas_r = 1.0;        oflSD_data.gasmed_r = 1.0;
    oflSD_data.lat_gas_r = 1.0;        oflSD_data.gaslat_r = 1.0;        oflSD_data.gaslat140_r = 1.0;
    oflSD_data.soleus_r  = 1.0; 
    oflSD_data.tib_post_r= 0.5;        oflSD_data.tibpost_r = 0.5;
    
    oflSD_data.flex_dig_r = 1.1;       oflSD_data.fdl_r = 1.1; 
    oflSD_data.flex_hal_r = 1.3;       oflSD_data.fhl_r = 1.3; 
    oflSD_data.tib_ant_r  = 0.8;       oflSD_data.tibant_r = 0.8; 
    oflSD_data.per_brev_r = 0.7;       oflSD_data.perbrev_r = 0.7; 
    oflSD_data.per_long_r = 0.6;       oflSD_data.perlong_r = 0.6; 
    oflSD_data.per_tert_r = NaN;                                    
    oflSD_data.ext_dig_r  = 1.1;        oflSD_data.edl_r = 1.1;
    oflSD_data.ext_hal_r  = 1.1;        oflSD_data.ehl_r = 1.1; 
    oflSD_data.ercspn_r  = NaN;
    oflSD_data.intobl_r  = NaN;
    oflSD_data.extobl_r  = NaN;    
    oflSD_data.add_long_r= 2.0;       oflSD_data.addlong_r = 2.0;
    oflSD_data.add_brev_r= 1.4;       oflSD_data.addbrev_r = 1.4;
    
    
    oflSD_data.glut_med1_l = 1.6;      oflSD_data.glmed1_l = 1.6;
    oflSD_data.glut_med2_l = 1.6;      oflSD_data.glmed2_l = 1.6;
    oflSD_data.glut_med3_l = 1.6;      oflSD_data.glmed3_l = 1.6;
    oflSD_data.glut_min1_l = NaN;      oflSD_data.glmin1_l = NaN;
    oflSD_data.glut_min2_l = NaN;      oflSD_data.glmin2_l = NaN;
    oflSD_data.glut_min3_l = NaN;      oflSD_data.glmin3_l = NaN;
    oflSD_data.semimem_l   = 1.8; 
    oflSD_data.semiten_l   = 4.1;   
    
    oflSD_data.bifemlh_l   = 2.6;      oflSD_data.bflh_l = 2.6;         oflSD_data.bflh140_l = 2.6; 
    oflSD_data.bifemsh_l   = 2.1;      oflSD_data.bfsh_l = 2.1;          oflSD_data.bfsh140_l = 2.1;  
    
    oflSD_data.sar_l       = 4.6;      oflSD_data.sart_l = 4.6;
    oflSD_data.add_mag1_l  = 3.4;      oflSD_data.addmagDist_l = 3.4;
    oflSD_data.add_mag2_l  = 3.0;      oflSD_data.addmagIsch_l = 3.0;
    oflSD_data.add_mag3_l  = 2.6;      oflSD_data.addmagMid_l  = 2.6;     
                                       oflSD_data.addmagProx_l = 2.0;
    oflSD_data.tfl_l       = NaN;
    oflSD_data.pect_l      = NaN;
    oflSD_data.grac_l      = 4.4;
    oflSD_data.glut_max1_l = 2.4;      oflSD_data.glmax1_l = 2.4;
    oflSD_data.glut_max2_l = 2.6;      oflSD_data.glmax2_l = 2.6;
    oflSD_data.glut_max3_l = 2.7;      oflSD_data.glmax3_l = 2.7;
    
    oflSD_data.iliacus_l   = 1.9;
    oflSD_data.psoas_l     = 1.7;
    oflSD_data.quad_fem_l  = NaN;
    oflSD_data.gem_l       = NaN;
    oflSD_data.peri_l      = NaN;      oflSD_data.piri_l = NaN;
    
    oflSD_data.rect_fem_l  = 1.3;      oflSD_data.recfem_l = 1.3;
    oflSD_data.vas_med_l = 2.3;        oflSD_data.vasmed_l = 2.3;
    oflSD_data.vas_int_l = 2.0;        oflSD_data.vasint_l = 2.0;
    oflSD_data.vas_lat_l = 1.8;        oflSD_data.vaslat_l = 1.8;        oflSD_data.vaslat140_l = 1.8; 
    
    oflSD_data.med_gas_l = 1.0;        oflSD_data.gasmed_l = 1.0;
    oflSD_data.lat_gas_l = 1.0;        oflSD_data.gaslat_l = 1.0;        oflSD_data.gaslat140_l = 1.0;
    oflSD_data.soleus_l  = 1.0; 
    oflSD_data.tib_post_l= 0.5;        oflSD_data.tibpost_l = 0.5;
    
    oflSD_data.flex_dig_l = 1.1;       oflSD_data.fdl_l = 1.1; 
    oflSD_data.flex_hal_l = 1.3;       oflSD_data.fhl_l = 1.3; 
    oflSD_data.tib_ant_l  = 0.8;       oflSD_data.tibant_l = 0.8; 
    oflSD_data.per_brev_l = 0.7;       oflSD_data.perbrev_l = 0.7; 
    oflSD_data.per_long_l = 0.6;       oflSD_data.perlong_l = 0.6; 
    oflSD_data.per_tert_l = NaN;                                    
    oflSD_data.ext_dig_l  = 1.1;        oflSD_data.edl_l = 1.1;
    oflSD_data.ext_hal_l  = 1.1;        oflSD_data.ehl_l = 1.1; 
    oflSD_data.ercspn_l  = NaN;
    oflSD_data.intobl_l  = NaN;
    oflSD_data.extobl_l  = NaN;    
    oflSD_data.add_long_l= 2.0;       oflSD_data.addlong_l = 2.0;
    oflSD_data.add_brev_l= 1.4;       oflSD_data.addbrev_l = 1.4;

    
    data = zeros(length(muscleNames),1);
    for i = 1:length(muscleNames)
         if isfield(oflSD_data,muscleNames{i})
            data(i,1) = oflSD_data.(muscleNames{i});
         else
             data(i,1) = NaN;
             disp(['Information about ' muscleNames{i}...
                 '  not found in script. Assigned value NaN']);
         end
    end    
end   