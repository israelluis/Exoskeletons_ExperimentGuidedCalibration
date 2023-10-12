% This function returns the values of the optimal fibers
% The data come from Rajagopal et al. (2016). (Ward et al.)
% We used 0 when no data were available.
%
% Author: Israel Luis
% Date: 04/01/2023
%     
function data = getOptimalFiber(muscleNames)

    % right side
    ofl_data.glut_med1_r = 7.3;      ofl_data.glmed1_r = 7.3;
    ofl_data.glut_med2_r = 7.3;      ofl_data.glmed2_r = 7.3;
    ofl_data.glut_med3_r = 7.3;      ofl_data.glmed3_r = 7.3;
    ofl_data.glut_min1_r = 6.8;      ofl_data.glmin1_r = 6.8; % from rajagopal
    ofl_data.glut_min2_r = 5.6;      ofl_data.glmin2_r = 5.6; % from rajagopal
    ofl_data.glut_min3_r = 3.8;      ofl_data.glmin3_r = 3.8; % from rajagopal
    ofl_data.semimem_r   = 6.9; %
    ofl_data.semiten_r   =19.3; %  
    
    ofl_data.bifemlh_r   = 9.8;      ofl_data.bflh_r = 9.8;          ofl_data.bflh140_r = 9.8; 
    ofl_data.bifemsh_r   =11.0;      ofl_data.bfsh_r =11.0;          ofl_data.bfsh140_r =11.0;  
    
    ofl_data.sar_r       =40.3;      ofl_data.sart_r =40.3;
    ofl_data.add_mag1_r  =17.7;      ofl_data.addmagDist_r =17.7; % from rajagopal
    ofl_data.add_mag2_r  =15.6;      ofl_data.addmagIsch_r =15.6; % from rajagopal
    ofl_data.add_mag3_r  =13.8;      ofl_data.addmagMid_r  =13.8; % from rajagopal     
                                     ofl_data.addmagProx_r =10.6; % from rajagopal
    ofl_data.tfl_r       = 9.5; % from rajagopal
    ofl_data.pect_r      = 0.0;
    ofl_data.grac_r      =22.8;
    ofl_data.glut_max1_r =14.7;      ofl_data.glmax1_r =14.7; % from rajagopal
    ofl_data.glut_max2_r =15.7;      ofl_data.glmax2_r =15.7; % from rajagopal
    ofl_data.glut_max3_r =16.7;      ofl_data.glmax3_r =16.7; % from rajagopal
    
    ofl_data.iliacus_r   =10.7;
    ofl_data.psoas_r     =11.7;
    ofl_data.quad_fem_r  = 0.0;
    ofl_data.gem_r       = 0.0;
    ofl_data.peri_r      = 2.6;      ofl_data.piri_r = 2.6; % from rajagopal
    
    ofl_data.rect_fem_r= 7.6;        ofl_data.recfem_r = 7.6;
    ofl_data.vas_med_r = 9.7;        ofl_data.vasmed_r = 9.7;
    ofl_data.vas_int_r = 9.9;        ofl_data.vasint_r = 9.9;
    ofl_data.vas_lat_r = 9.9;        ofl_data.vaslat_r = 9.9;        ofl_data.vaslat140_r = 9.9; 
    
    ofl_data.med_gas_r = 5.1;        ofl_data.gasmed_r = 5.1;
    ofl_data.lat_gas_r = 5.8;        ofl_data.gaslat_r = 5.8;        ofl_data.gaslat140_r = 5.8;
    ofl_data.soleus_r  = 4.4; 
    ofl_data.tib_post_r= 3.8;        ofl_data.tibpost_r= 3.8;
    
    ofl_data.flex_dig_r = 4.5;       ofl_data.fdl_r = 4.5; 
    ofl_data.flex_hal_r = 5.3;       ofl_data.fhl_r = 5.3; 
    ofl_data.tib_ant_r  = 6.8;       ofl_data.tibant_r = 6.8; % 
    ofl_data.per_brev_r = 4.5;       ofl_data.perbrev_r = 4.5; 
    ofl_data.per_long_r = 5.1;       ofl_data.perlong_r = 5.1; 
    ofl_data.per_tert_r = 0.0;                                    
    ofl_data.ext_dig_r  = 6.9;       ofl_data.edl_r = 6.9;
    ofl_data.ext_hal_r  = 7.5;       ofl_data.ehl_r = 7.5; 
    ofl_data.ercspn_r  = 0.0;
    ofl_data.intobl_r  = 0.0;
    ofl_data.extobl_r  = 0.0;    
    ofl_data.add_long_r=10.8;        ofl_data.addlong_r =10.8;
    ofl_data.add_brev_r=10.3;        ofl_data.addbrev_r =10.3;
    
    % left side
    ofl_data.glut_med1_l = 7.3;      ofl_data.glmed1_l = 7.3;
    ofl_data.glut_med2_l = 7.3;      ofl_data.glmed2_l = 7.3;
    ofl_data.glut_med3_l = 7.3;      ofl_data.glmed3_l = 7.3;
    ofl_data.glut_min1_l = 6.8;      ofl_data.glmin1_l = 6.8; % from rajagopal
    ofl_data.glut_min2_l = 5.6;      ofl_data.glmin2_l = 5.6; % from rajagopal
    ofl_data.glut_min3_l = 3.8;      ofl_data.glmin3_l = 3.8; % from rajagopal
    ofl_data.semimem_l   = 6.9; %
    ofl_data.semiten_l   =19.3; %  
    
    ofl_data.bifemlh_l   = 9.8;      ofl_data.bflh_l = 9.8;          ofl_data.bflh140_l = 9.8; 
    ofl_data.bifemsh_l   =11.0;      ofl_data.bfsh_l =11.0;          ofl_data.bfsh140_l =11.0;  
    
    ofl_data.sar_l       =40.3;      ofl_data.sart_l =40.3;
    ofl_data.add_mag1_l  =17.7;      ofl_data.addmagDist_l =17.7; % from rajagopal
    ofl_data.add_mag2_l  =15.6;      ofl_data.addmagIsch_l =15.6; % from rajagopal
    ofl_data.add_mag3_l  =13.8;      ofl_data.addmagMid_l  =13.8; % from rajagopal     
                                     ofl_data.addmagProx_l =10.6; % from rajagopal
    ofl_data.tfl_l       = 9.5; % from rajagopal
    ofl_data.pect_l      = 0.0;
    ofl_data.grac_l      =22.8;
    ofl_data.glut_max1_l =14.7;      ofl_data.glmax1_l =14.7; % from rajagopal
    ofl_data.glut_max2_l =15.7;      ofl_data.glmax2_l =15.7; % from rajagopal
    ofl_data.glut_max3_l =16.7;      ofl_data.glmax3_l =16.7; % from rajagopal
    
    ofl_data.iliacus_l   =10.7;
    ofl_data.psoas_l     =11.7;
    ofl_data.quad_fem_l  = 0.0;
    ofl_data.gem_l       = 0.0;
    ofl_data.peri_l      = 2.6;      ofl_data.piri_l = 2.6; % from rajagopal
    
    ofl_data.rect_fem_l= 7.6;        ofl_data.recfem_l = 7.6;
    ofl_data.vas_med_l = 9.7;        ofl_data.vasmed_l = 9.7;
    ofl_data.vas_int_l = 9.9;        ofl_data.vasint_l = 9.9;
    ofl_data.vas_lat_l = 9.9;        ofl_data.vaslat_l = 9.9;        ofl_data.vaslat140_l = 9.9; 
    
    ofl_data.med_gas_l = 5.1;        ofl_data.gasmed_l = 5.1;
    ofl_data.lat_gas_l = 5.8;        ofl_data.gaslat_l = 5.8;        ofl_data.gaslat140_l = 5.8;
    ofl_data.soleus_l  = 4.4; 
    ofl_data.tib_post_l= 3.8;        ofl_data.tibpost_l= 3.8;
    
    ofl_data.flex_dig_l = 4.5;       ofl_data.fdl_l = 4.5; 
    ofl_data.flex_hal_l = 5.3;       ofl_data.fhl_l = 5.3; 
    ofl_data.tib_ant_l  = 6.8;       ofl_data.tibant_l = 6.8; % 
    ofl_data.per_brev_l = 4.5;       ofl_data.perbrev_l = 4.5; 
    ofl_data.per_long_l = 5.1;       ofl_data.perlong_l = 5.1; 
    ofl_data.per_tert_l = 0.0;                                    
    ofl_data.ext_dig_l  = 6.9;       ofl_data.edl_l = 6.9;
    ofl_data.ext_hal_l  = 7.5;       ofl_data.ehl_l = 7.5; 
    ofl_data.ercspn_l  = 0.0;
    ofl_data.intobl_l  = 0.0;
    ofl_data.extobl_l  = 0.0;    
    ofl_data.add_long_l=10.8;        ofl_data.addlong_l =10.8;
    ofl_data.add_brev_l=10.3;        ofl_data.addbrev_l =10.3;
    
    data = zeros(length(muscleNames),1);
    for i = 1:length(muscleNames)
         if isfield(ofl_data,muscleNames{i})
            data(i,1) = ofl_data.(muscleNames{i});
         else
             data(i,1) = 0.0;
             disp(['Information about ' muscleNames{i}...
                 '  not found in script. Assigned value 0']);
         end
    end    
end   