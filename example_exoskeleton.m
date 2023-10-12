addpath(genpath('C:\Users\israe\Dropbox\PhDWork\Simulation\code\MuscleRedundancySolver-master\Casadi_installed'));
addpath(genpath('C:\Users\israe\Dropbox\PhDWork\Simulation\code\Exo&Calibration'));
clc; close all; clear Misc
%% Input information
currentFolder = pwd;

% paths and folders
ExampleFolder= 'DataExample';
DataFolder   = 'DataDigitalized';
OutPath      = fullfile(currentFolder,'Results_Exo');  

% experimental data
Misc.IKfile = {fullfile(currentFolder,ExampleFolder,'IK.mot')};
Misc.IDfile = {fullfile(currentFolder,ExampleFolder,'ID.sto')};
model_path  =  fullfile(currentFolder,ExampleFolder,'Model.osim');

% select the leg's side
Misc.side_sel      = 'r';

% input dofs - select the DOFs you want to include in the optimization
Misc.DofNames_Input={['hip_flexion_' Misc.side_sel] ['hip_adduction_' Misc.side_sel] ['hip_rotation_' Misc.side_sel] ['knee_angle_' Misc.side_sel]  ['ankle_angle_' Misc.side_sel]}; 

% ranges
time = [0 10];

% Time of toe_off helps set the limits of the timing for exoskeleton
% asssistance
Misc.toe_off_timing= 1.31;

% Original FMo values from Rajagopal might be to high to see a muscle
% activation reduction. I decreased it here
Misc.updParams.FMo.names      ={'soleus' 'gasmed' 'gaslat'};
Misc.updParams.FMo.values     =[3549 1558 683]; % i.e., this are the new values

% label of muscle-tendon parameters
param_label={'FMo' 'lMo' 'lTs' 'alphao' 'vMmax' 'kpe' 'so' 'sM'};
Misc.param_label   =param_label;

% run muscle analysis
Misc.GetAnalysis = 1;

% We will use previously optimized lMo, lTs, kT, and passive force-length
% curves. Run "" to have access to the optimized values
Misc.Param_read    = {'selected'};
Misc.myParams=labelParams(Results_all.kT.calibratedMRS,Results_all.params.calibratedMRS,param_label);

% Weight terms in the objective function (otherwise it uses default values)
Misc.wEMG   =    0.00; 
Misc.wAct   =    1.00; 
Misc.wTres  = 1000.00; 
Misc.wVm    =    0.001; 

% WORKFLOW SETUP
% Here we setup 4 workflows. We run the muscle redundacy with assistive
% moments
% normal: no assistive moment
% Ankle_spring: Assistive moment with clutchable spring to support plantarflexion
% Ankle_motor: Assistive moment with ideal actuation to support plantarflexion
% Knee_custom: Insert your own assistive moment to support knee extension
% This last feature isnt discussed in the paper, but it is just fun to have
% it for quickly exploring alternative research questions
Misc.muscleFiberCal= 0; 

workFlow_Normal       =1;
workFlow_Ankle_spring =1;
workFlow_Ankle_motor  =1;
workFlow_Knee_custom  =1;

Misc.Exo_Enable    = 0; % no exo
if workFlow_Normal==1
Misc.OutName= 'normal';

[Results_normal,DatStore,Misc] = solveMuscleRedundancy_ExoCal(model_path,time,OutPath,Misc);
end

Misc.Exo_Enable    = 1; % exo enabled
if workFlow_Ankle_spring==1
Misc.Exo_Mode="optimal"; Misc.Exo_Type="quasi"; % optimal assistance with spring
Misc.Exo_Dof={['ankle_angle_' Misc.side_sel]}; Misc.Exo_group=-1; % ankle plantarflexion - / ankle dorsiflexion +
Misc.OutName= 'ankle_spring';

[Results_ankSpring,DatStore,Misc] = solveMuscleRedundancy_ExoCal(model_path,time,OutPath,Misc);
end

if workFlow_Ankle_spring==1
Misc.Exo_Mode="optimal"; Misc.Exo_Type="quasi"; % optimal assistance with spring
Misc.Exo_Dof={['ankle_angle_' Misc.side_sel]}; Misc.Exo_group=-1; % ankle plantarflexion - / ankle dorsiflexion +
Misc.OutName= 'ankle_spring';

[Results_ankSpring,DatStore,Misc] = solveMuscleRedundancy_ExoCal(model_path,time,OutPath,Misc);
end

if workFlow_Ankle_motor==1
Misc.Exo_Mode='optimal'; Misc.Exo_Type="active"; % optimal assistance with motor
Misc.Exo_Dof={['ankle_angle_' Misc.side_sel]}; Misc.Exo_group=-1; % ankle plantarflexion - / ankle dorsiflexion +
Misc.OutName= 'ankle_motor';

[Results_ankMotor,DatStore,Misc] = solveMuscleRedundancy_ExoCal(model_path,time,OutPath,Misc);
end

if workFlow_Knee_custom==1
Misc.Exo_Mode='optimal'; Misc.Exo_Type="active"; % optimal assistance with motor
Misc.Exo_Dof={['knee_angle_' Misc.side_sel]}; Misc.Exo_group=-1; % knee extension - / knee flexion +
Misc.OutName= 'knee_custom';

% generate a torque profile
% any moment can be added to Misc.torqueProfile. For simplicity I use a sin function
t_exo=[5 20]; % percentage of gait cycle you want to create the torque
t_d=t_exo(end)-t_exo(1)+1;
ang=linspace(0,pi,t_d);
torqueSine=sin(ang);

peak= -50;
torqueAssistive=zeros(1,100);
torqueAssistive(t_exo(1):t_exo(2))=torqueSine.*peak;
Misc.torqueProfile(1)={torqueAssistive};

[Results_kneCustom,DatStore,Misc] = solveMuscleRedundancy_ExoCal(model_path,time,OutPath,Misc);
end
%% Ploting all results with assistive moments
state_list={'MActivation' 'TForce'};
yaxis_values ={[0 1] [0 1500]};
title_list={'muscle activations in normal and assisted conditions' 'muscle-tendon forces in normal and assisted conditions'};
yaxis_list={'muscle activations ' 'muscle-tendon forces [N]'};

MuscleNames= DatStore.MuscleNames;
NMuscle    = length(MuscleNames);

for state_opt=1:length(state_list)
fig=figure(state_opt); set(gcf,'color','w','Visible','on'); clf;   

for mus_sel=1:NMuscle
    muscle_name=MuscleNames{mus_sel};
    normal =Results_normal.(state_list{state_opt}).genericMRS(mus_sel,1+extra_frames:end-extra_frames);
    case_1 =Results_ankSpring.(state_list{state_opt}).genericMRS(mus_sel,1+extra_frames:end-extra_frames);
    case_2 =Results_ankMotor.(state_list{state_opt}).genericMRS(mus_sel,1+extra_frames:end-extra_frames);
    case_3 =Results_kneCustom.(state_list{state_opt}).genericMRS(mus_sel,1+extra_frames:end-extra_frames);
    
    x_axis      =linspace(0,100,length(normal));
    
    subplot(5,8,mus_sel)
    hold on;
    plot(x_axis,normal,'color',	"#4c4c4c",'lineWidth',3,'lineStyle','-');
    plot(x_axis,case_1,'color',	"#0072BD",'lineWidth',3,'lineStyle','-');
    plot(x_axis,case_2,'color',	"#A2142F",'lineWidth',3,'lineStyle','-');
    plot(x_axis,case_3,'color', "#EDB120",'lineWidth',3,'lineStyle','-.');
    ylim(yaxis_values{state_opt});
    set(gca,'FontSize',10);
    title(muscle_name)
end
sgtitle(title_list(state_opt),'fontSize',20)

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
label_y = ylabel(han,yaxis_list{state_opt},'FontSize',17); label_y.Position(1) = -0.05; label_y.Position(2) = 0.5;
label_x = xlabel(han,'gait cycle [%]'       ,'FontSize',17);   label_x.Position(1) = 0.5; label_x.Position(2) = -0.05;
end
