addpath(genpath('C:\Users\israe\Dropbox\PhDWork\Simulation\code\MuscleRedundancySolver-master\Casadi_installed'));
addpath(genpath('C:\Users\israe\Dropbox\PhDWork\Simulation\code\Exo&Calibration'));
clc; close all; clear Misc
% Input information
currentFolder = pwd;

% paths and folders
ExamplePath= 'DataExample';
DataPath   = 'DataDigitalized';
OutPath    = fullfile(currentFolder,'Results_PassMom');  

% model:rajagopal: knee bending positive, 2392: knee bending negative angle
model_conf ={'rajagopal'};               Misc.model_conf=model_conf; 

% select the leg's side
Misc.side_sel      = 'r';

% experimental data
ank_file =fullfile(DataPath,['TrialAnk_' Misc.side_sel '_kne0hip0_kneBen15hip0_kneBen60hip0_' model_conf{:} 'Model']);
hip_file =fullfile(DataPath,['TrialHip_' Misc.side_sel '_kneBen15ank0_kneBen60ank0_' model_conf{:} 'Model']);
kne_file =fullfile(DataPath,['TrialKne_' Misc.side_sel '_ankDor20hip0_ankPla15hip0_ankDor20hipExt15_' model_conf{:} 'Model']);
neu_file =fullfile(DataPath,['TrialNeu_' Misc.side_sel '_hipExt7hipAbd2kne0ankPla40_allModels']);

Misc.IKfile = {[ank_file '.mot'] [hip_file '.mot'] [kne_file '.mot'] [neu_file '.mot']};
Misc.IDfile = {[ank_file '.sto'] [hip_file '.sto'] [kne_file '.sto'] [neu_file '.sto']};

model_path  =  fullfile(currentFolder,ExamplePath,'Model.osim');

Misc.model_path = model_path;
Misc.OutPath    = OutPath;

% input dofs - select the DOFs you want to include in the optimization
Misc.DofNames_Input={{['ankle_angle_' Misc.side_sel]};{['hip_flexion_' Misc.side_sel]};{['knee_angle_' Misc.side_sel]};...        % for passive forces
{['hip_flexion_' Misc.side_sel] ['hip_adduction_' Misc.side_sel] ['hip_rotation_' Misc.side_sel] ['knee_angle_' Misc.side_sel] ['ankle_angle_' Misc.side_sel]}};   % for neutral

% ranges
time = [0 10; 0 10; 0 10; 0 10];

% label of muscle-tendon parameters
param_label={'FMo' 'lMo' 'lTs' 'alphao' 'vMmax' 'kpe' 'so' 'sM'};
Misc.param_label   =param_label;

% run muscle analysis
Misc.RunAnalysis = 1;

% Adjust semimem
% semimen's tendon needs to be adjusted as fiber length is, otherwise, too short. 
% You can decrease the tendon slack length (a) based on a small percentage or
% (b) recompute it to be non-negative while muscle-tendon unit is shortening.

Misc.semimem_adjustment=2; % change to "2" to have the lTs recomputed
if Misc.semimem_adjustment==1 
    Misc.updParams.lTs.names      ={'semimem'};
    Misc.updParams.lTs.percentages=[0.98]; % i.e., 0.98 its original value
end

% WORKFLOW SETUP 
% Here we setup 2 workflows.
% Generic:     Generic parameters, no recorded passive angle-moment
% Calibration: Calibration of parameters using recorded passive angle-moment

workFlow_generic     =1;
workFlow_calibration =1;

if workFlow_generic==1
    Misc.Param_read={'model'};
    Misc.Mode_opti={'none'}; % no optimization
    Misc.OutName = ['GenericPass_Results'];
    
    Misc.add_Passconstraints=0;
    [Results_PassGen,DatStore_PassGen,Misc_PassGen] = solveMuscleRedundancy_calibrationPassive([],time,[],Misc);
end

if workFlow_calibration==1
    Misc.Param_read={'model'};
    Misc.Mode_opti={'passive'}; % passive force-length curve optimization
    Misc.OutName = ['GenericPass_Results'];
    
    Misc.add_Passconstraints=1;
    [Results_PassCal,DatStore_PassCal,Misc_PassCal] = solveMuscleRedundancy_calibrationPassive([],time,[],Misc);
end
%% PLOT PASSIVE FORCES
nTrial_per_file =[3 2 3];
joint_order     =[1 3 2];
jointLeg_label   = {'ankle','knee','hip'};
description_label= {'knee flexion  0°' 'knee flexion 15°' 'knee flexion 60°';
                   {'ankle dorsiflexion 20°'; 'hip neutral 0°'} {'ankle plantarflexion 15°'; 'hip neutral 0°'} {'ankle dorsiflexion 20°'; 'hip extension 15°'};
                   'knee flexion 15°' 'knee flexion 60°' ' '};
jointAngle_label = {'ankle dorsiflexion [deg]' 'knee flexion [deg]' 'hip extension [deg]'};

lim_x_list={[-30 20]; [0 75]; [-15 40]};
lim_y_list={[-45 45]; [-55 35]; [-50 40]}; %[-45 20]; [-60 30]; [-50 40]
color_pass={'#191919' '#4682B4' '#cc0000'};
frame_per_trial=50;   
fig=figure(3); clf;  set(gcf,'color','w','Visible','on','Position',[50 50 1000 700]);
for joint_opt=1:3
    joint_sel=joint_order(joint_opt);
    for trial_sel=1:nTrial_per_file(joint_sel)
        
        sel_off=(trial_sel-1)*frame_per_trial;
        q_exp_array= DatStore_PassCal(joint_sel).q_exp(1+sel_off:frame_per_trial+sel_off);
        T_exp_array= DatStore_PassCal(joint_sel).T_exp(1+sel_off:frame_per_trial+sel_off);
        T_gen_array= Results_PassGen.MMoment(joint_sel).genericMRS(1+sel_off:frame_per_trial+sel_off);
        T_cal_array= Results_PassCal.MMoment(joint_sel).calibratedMRS(1+sel_off:frame_per_trial+sel_off);
        
        rcor_gen = corr(T_exp_array, T_gen_array);
        rmse_gen = sqrt(mean((T_exp_array - T_gen_array).^2));
       
        rcor_cal = corr(T_exp_array, T_cal_array);
        rmse_cal = sqrt(mean((T_exp_array - T_cal_array).^2));
       
        subplot(3,3,joint_opt+(trial_sel-1)*3)
        hold on;
        plot(q_exp_array,T_exp_array,'o','MarkerEdgeColor',"#000000",'MarkerFaceColor',color_pass{1},'markerSize',5);
        plot(q_exp_array,T_gen_array,'o','MarkerEdgeColor',"#000000",'MarkerFaceColor',color_pass{2},'markerSize',5);
        plot(q_exp_array,T_cal_array,'o','MarkerEdgeColor',"#000000",'MarkerFaceColor',color_pass{3},'markerSize',5);
       
        text(lim_x_list{joint_opt}(1)+7,lim_y_list{joint_opt}(2)-2,['GEN r=' num2str(rcor_gen,'%4.2f') ' RMSE=' num2str(rmse_gen,'%4.2f')],'fontSize',8,'color',color_pass{2},'FontWeight','bold','HorizontalAlignment','left')
        text(lim_x_list{joint_opt}(1)+7,lim_y_list{joint_opt}(2)-12,['CAL r=' num2str(rcor_cal,'%4.2f') ' RMSE=' num2str(rmse_cal,'%4.2f')],'fontSize',8,'color',color_pass{3},'FontWeight','bold','HorizontalAlignment','left')
       
        text(lim_x_list{joint_opt}(1)+2,lim_y_list{joint_opt}(1)+20,description_label{joint_opt+(trial_sel-1)*3},'fontSize',10)
        xlim(lim_x_list{joint_opt,:});
        ylim(lim_y_list{joint_opt,:});
       
        set(gca,'FontSize',12);
        xlabel(jointAngle_label{joint_opt},'FontSize',12);
        ylabel('moment [Nm]','FontSize',12)
       
       
        if trial_sel==1 && (joint_opt==1 || joint_opt==2 || joint_opt==3)
            title(jointLeg_label{joint_opt},'FontSize',20);
        end
    end
end

subplot(3,3,9)
hold on;
plot(0,0,'o','MarkerEdgeColor',"#000000",'MarkerFaceColor',color_pass{1},'markerSize',10);
plot(0,0,'o','MarkerEdgeColor',"#000000",'MarkerFaceColor',color_pass{2},'markerSize',10); hold on;
plot(0,0,'o','MarkerEdgeColor',"#000000",'MarkerFaceColor',color_pass{3},'markerSize',10); hold on;

axis([10,11,10,11]) %move dummy points out of view
legend('experimental','generic','calibrated','location','southwest','fontSize',15); legend boxoff  
axis off %hide axis