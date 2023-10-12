function [Results,DatStore,Misc] = solveMuscleRedundancy_calibrationPassive(model_path,time,OutPath,Misc)
%% Run muscle tendon estimator:
% No filter of data (as we assume that velocity = 0 )
Misc.f_order_dM = NaN;
Misc.f_order_lMT= NaN;
Misc.f_order_IK = NaN;
Misc.f_order_ID = NaN;

% update default settings
Misc = DefaultSettings_Gupta(Misc);

% read the muscle properties
[Misc] = getMuscleProperties(Misc.model_path,Misc);

% number of motion trials
Misc.nTrials = length(Misc.IKfile);

%check if we have to adapt the start and end time so that it corresponds to
%time frames in the IK solution
[time] = Check_TimeIndices(Misc,time);
Misc.time=time;

%% Perform muscle analysis for all trials
DatStore = struct;
MuscleAnalysisPath=fullfile(Misc.OutPath,'MuscleAnalysis');
Misc.MuscleAnalysisPath=MuscleAnalysisPath;
if ~exist(MuscleAnalysisPath,'dir')
    mkdir(MuscleAnalysisPath);
end
for i = 1:Misc.nTrials
    % select the IK and ID file
    IK_path_trial = Misc.IKfile{i};
    % Run muscle analysis
    if Misc.RunAnalysis
        disp('MuscleAnalysis Running .....');
        OpenSim_Muscle_Analysis(IK_path_trial,Misc.model_path,MuscleAnalysisPath,[time(i,1) time(i,end)],Misc.DofNames_Input{i,1})
        disp('MuscleAnalysis Finished');
    end
end

%% Extract muscle information
% Get number of degrees of freedom (dofs), muscle-tendon lengths, moment
% arms, stiffness and shift for the selected muscles.
for trial = 1:Misc.nTrials
    [~,Misc.MAtrialName{trial},~]=fileparts(Misc.IKfile{trial});
end

% select muscles with a moment arm for the input dofs
Misc = getMuscles4DOFS_Gupta(Misc);

% get IK, ID, muscle analysis data
[Misc,DatStore] = getMuscleInfo_GuptaNoFil(Misc,DatStore);

% read the muscle properties
[Misc.params,Misc.lMo,Misc.lTs,Misc.FMo,Misc.alphao,Misc.kT,Misc.shift]=ReadMuscleParameters_upd(Misc.model_path,DatStore(4).MuscleNames,Misc);
% complete for both legs
[row,col]=size(Misc.params); totalMuscles=length(DatStore(4).MuscleNames);
if Misc.side_sel=='l'
    if col == totalMuscles
    else    
        Misc.params(:,1:80)=[zeros(8,40) Misc.params];
        Misc.kT(:,1:80)    =[zeros(1,40) Misc.kT];
        Misc.shift(:,1:80) =[zeros(1,40) Misc.shift];
    end
elseif Misc.side_sel=='r'
    if col == totalMuscles
    else
        Misc.params(:,1:80)=[Misc.params zeros(8,40)];
        Misc.kT(:,1:80)    =[Misc.kT     zeros(1,40)];
        Misc.shift(:,1:80) =[Misc.shift  zeros(1,40)];   
    end
end

% display warnings in muscle selection
[Misc] = Warnings_MuscleNames_all(DatStore,Misc);

% get indexes of the muscles for which optimal fiber length, tendon stiffness are estimated
% [DatStore,Misc] = GetIndices(DatStore,Misc);

% get the EMG and ultrasound information
% [Misc,DatStore] = GetEMGInfo(Misc,DatStore);
% [Misc,DatStore] = GetUSInfo(Misc,DatStore);

% EMG_data = ReadMotFile(Misc.EMGfile{1});
% US_data  = ReadMotFile(Misc.USfile{1});
    
% get the number of muscles in a vector
NMuscles = zeros(Misc.nTrials,1);
for trial = 1:Misc.nTrials
    NMuscles(trial) = DatStore(trial).nMuscles;
end

%% set solver
% Create an NLP solver
SolverSetup.nlp.solver = 'ipopt';
SolverSetup.derivatives.derivativelevel = 'second';
SolverSetup.optionssol.ipopt.nlp_scaling_method = 'gradient-based';
SolverSetup.optionssol.ipopt.linear_solver = 'mumps';
SolverSetup.optionssol.ipopt.tol = 1e-6;
SolverSetup.optionssol.ipopt.max_iter = 10000;
%% initial guess
% initial guess is based on rigid tendon (and ignores pennation angle)
plot_intermediate=0;
for trial = 1:Misc.nTrials
    nMuscles = length(DatStore(trial).MuscleNames);
    Nfr = length(DatStore(trial).time);
    lTs = Misc.params(3,Misc.idx_allMuscleList{trial})';
    IG(trial).lM_projected = zeros(nMuscles,Nfr);
    IG(trial).lMtilde      = zeros(nMuscles,Nfr);
    IG(trial).a            = zeros(nMuscles,Nfr);

    lMT     = DatStore(trial).LMT(:,:)';
    lMGuess = lMT-lTs;
    
    if Misc.semimem_adjustment==2
        for mus=1:nMuscles
            if ~isempty(find(lMGuess(mus,:)<0))
                lTs_upd=min(lMGuess(mus,:));
                lTs_new=lTs(mus)+lTs_upd-0.001;
                Misc.params(3,Misc.idx_allMuscleList{trial}(mus))=lTs_new;
                lMGuess(mus,:) = lMT(mus,:)-lTs_new;
            end
        end
    end
    
    lMo = Misc.params(2,Misc.idx_allMuscleList{trial})';
    alphao = Misc.params(4,Misc.idx_allMuscleList{trial})';
    w = lMo.*sin(alphao);
    IG(trial).lM_projected(:,:) = sqrt((lMGuess.^2 - w.^2));
    IG(trial).lMtilde(:,:) = lMGuess./lMo;

    if plot_intermediate==1
        figure;
        plot(IG(trial).lMtilde','lineStyle','-','lineWidth',2); hold on;
        xlabel('frames');
        ylabel('rigid tendon estimate normalized fiber length');
        legend(DatStore(trial).MuscleNames)
    end
end
%% formulate optimization problem
addpath(genpath('C:\Users\israe\Dropbox\PhDWork\Simulation\code\myFunctions'));
% we mimize difference between inverse dynamic moments and moments
% generated by muscles. We account for tendon compliance here and solve
% for equilibrium between tendon force and the force in the passive element
% in the muscle fiber

% get constants
trial = 4;
MuscProperties.params = Misc.params(:,Misc.idx_allMuscleList{trial})';
lMo_list              = MuscProperties.params(:,2);
lTs_list              = MuscProperties.params(:,3);
kpe_list              = MuscProperties.params(:,6);
so_list               = MuscProperties.params(:,7);
sM_list               = MuscProperties.params(:,8);

MuscleNames_total   = DatStore(trial).MuscleNames;
lMo_ward            = getOptimalFiber(MuscleNames_total)/100;   % [m]
nMuscles            = length(MuscleNames_total);

import casadi.*
opti    = casadi.Opti();    % create opti structure

clear TR
J  = 0;

w_trial  = [1.0 1.0 1.0 0.0]; % 1.0(x3) 1.0
w_regular=  0.10; % 0.01

lMo_var = MX(nMuscles,1);
lTs_var = MX(nMuscles,1);
kpe_var = MX(nMuscles,1);
so_var  = MX(nMuscles,1);
sM_var  = MX(nMuscles,1);

if strcmp(Misc.Mode_opti,'none')
    lMo_var = lMo_list;
    lTs_var = lTs_list;
    kpe_var = kpe_list;
    so_var  = so_list;
    sM_var  = sM_list;
elseif strcmp(Misc.Mode_opti,'passive')
    lMo_per_lim = 0.00; % 0.01 (before) 0.10 makes large differences in normalized fibers
    lTs_per_lim = 0.00; % 0.01 (before)

    kpe_initial= 4.0; kpe_bound=1.00; 
    so_initial = 1.0; so_bound =0.25;
    sM_initial = 0.6; sM_bound =0.20;
    
    lMo_lim =[(1-lMo_per_lim).*lMo_list (1+lMo_per_lim).*lMo_list];
    lTs_lim =[(1-lTs_per_lim).*lTs_list (1+lTs_per_lim).*lTs_list];
    kpe_lim =[ones(nMuscles,1).*kpe_initial-kpe_bound  ones(nMuscles,1).*kpe_initial+kpe_bound];
    so_lim  =[ones(nMuscles,1).*so_initial-so_bound    ones(nMuscles,1).*so_initial+so_bound];
    sM_lim  =[ones(nMuscles,1).*sM_initial-sM_bound    ones(nMuscles,1).*sM_initial+sM_bound];

    % define variables
    lMo_var = opti.variable(nMuscles,1);
    lTs_var = opti.variable(nMuscles,1);
    kpe_var = opti.variable(nMuscles,1);
    so_var  = opti.variable(nMuscles,1);
    sM_var  = opti.variable(nMuscles,1);

    % constraints
    opti.subject_to(lMo_lim(:,1) < lMo_var < lMo_lim(:,2));             
    opti.subject_to(lTs_lim(:,1) < lTs_var < lTs_lim(:,2));
    opti.subject_to(kpe_lim(:,1) < kpe_var < kpe_lim(:,2)); 
    opti.subject_to( so_lim(:,1) < so_var  <  so_lim(:,2)); 
    opti.subject_to( sM_lim(:,1) < sM_var  <  sM_lim(:,2)); 

    % set initial guess
    opti.set_initial(lMo_var, lMo_list);
    opti.set_initial(lTs_var, lTs_list);
    opti.set_initial(kpe_var, kpe_initial);
    opti.set_initial(so_var , so_initial);
    opti.set_initial(sM_var , sM_initial);

    % regularization
    J= J + w_regular.*(sumsqr((kpe_var-kpe_initial)/100)+sumsqr((so_var -so_initial)/100)+...
                       sumsqr((sM_var -sM_initial)/100))/nMuscles;
end

MuscProperties_params_all =[MuscProperties.params(:,1) lMo_var lTs_var MuscProperties.params(:,4:5) kpe_var so_var sM_var];
%%
for trial=1:Misc.nTrials  
    % get number of muscles and frames
    nMuscles   = length(DatStore(trial).MuscleNames);
    Nfr        = length(DatStore(trial).time);
    selectedMus= DatStore(trial).MuscleNames;
    
    inds=zeros(nMuscles,1);
    for i=1:nMuscles
        inds(i) = find(strcmp(MuscleNames_total,selectedMus(i)));
    end
    
    % select muscle properties
    MuscProperties.params = Misc.params(:,Misc.idx_allMuscleList{trial})';
    MuscProperties.kT     = Misc.kT(:,Misc.idx_allMuscleList{trial}')';
    MuscProperties.shift  = Misc.shift(:,Misc.idx_allMuscleList{trial}')';
    alphao                = MuscProperties.params(:,4);

    % get variables
    lMo=lMo_var(inds,1);
    lTs=lTs_var(inds,1);
    kpe=kpe_var(inds,1);
    so = so_var(inds,1);
    sM = sM_var(inds,1);

    % muscle activation, bounds are given later
    TR(trial).a=zeros(nMuscles,Nfr);
    a=TR(trial).a;

    % updated for optimal fiber length and tendon slack length estimation
    MuscProperties_params =[MuscProperties.params(:,1:1) lMo lTs MuscProperties.params(:,4:5) kpe so sM];
        
    % projected fiber length as optimization variable
    TR(trial).lM_projected = opti.variable(nMuscles,Nfr);
    lM_projected=TR(trial).lM_projected;
    opti.subject_to(1e-4 < lM_projected(:)); % We impose that projected muscle fiber length has strict positive length
    opti.set_initial(lM_projected,IG(trial).lM_projected);

    % create reserve actuators
    TR(trial).aTk = opti.variable(DatStore(trial).nDOF,Nfr);
    aTk=TR(trial).aTk;
    opti.subject_to(-1< aTk < 1)

    % lMtilde as optimization variable, bounds are given later
    TR(trial).lMtilde = opti.variable(nMuscles,Nfr);
    lMtilde=TR(trial).lMtilde;
    opti.subject_to(0.05 < lMtilde < 2.00);
    opti.set_initial(lMtilde,IG(trial).lMtilde);
            
    % constraint on projected fiber length
    lM = lMtilde.*lMo;
    w = lMo.*sin(alphao);
    opti.subject_to(lM.^2 - w.^2 == lM_projected.^2);

    % find passive
    TR(trial).TsimVect = MX(Nfr,DatStore(trial).nDOF);
    
    for k=1:Nfr
        % solve muscle equilibrium
        [err, FTk] = ForceEquilibrium_lMtildeState_optPassive(a(:,k),lMtilde(:,k),...
            0,lM_projected(:,k),DatStore(trial).LMT(k,:)',MuscProperties_params,...
            MuscProperties.kT,MuscProperties.shift);

        % Moment constraints
        for dof = 1:DatStore(trial).nDOF
            T_exp = DatStore(trial).T_exp(k,dof); % + 1 due to the time col
            T_sim = FTk'*squeeze(DatStore(trial).dM(k,dof,:)) + Misc.Topt*aTk(dof,k);
            opti.subject_to(T_exp - T_sim == 0);
            TR(trial).TsimVect(k,dof) = FTk'*squeeze(DatStore(trial).dM(k,dof,:));
        end
        % impose equilibrium
        opti.subject_to(err == 0);
    end
    J = J + w_trial(trial)*(1000*sumsqr(aTk)/Nfr/DatStore(trial).nDOF+sumsqr(a)/Nfr/nMuscles);
end
   
opti.minimize(J);
opti.solver(SolverSetup.nlp.solver,SolverSetup.optionssol);

% begin - copy:constraints for calibrating parameters in bi-articular muscles
% end -> copy
additional_constraints=Misc.add_Passconstraints;
if additional_constraints==1
%% Additional constraints 
% muscle groups have the same values (since we for example cannot descriminate between parameters of the vastus
% lateral, medialis and intermdius we impose that the parameters should be the same)
if strcmp(Misc.Mode_opti,'passive')
muscle_coupling_passive=cell(1,Misc.nTrials);
for trial=1:Misc.nTrials
    for i=1:length(Misc.DofNames_Input{trial})
        dof_name=Misc.DofNames_Input{trial}{i}(1:end-2);%degree of freedom without side
        side    =Misc.DofNames_Input{trial}{i}(end);    %side only
        [muscle_coupling_passive{trial}]= muscleCouplingPassive(Misc.model_conf,dof_name,side);
    end
end

for trial=1:Misc.nTrials 
    [row,col]=size(muscle_coupling_passive{trial});
%     MuscleNames=DatStore(trial).MuscleNames;
    for i=1:row
        index1=find(contains(MuscleNames_total,muscle_coupling_passive{trial}(i,1)));
        index2=find(contains(MuscleNames_total,muscle_coupling_passive{trial}(i,2)));
        if ~isempty(index1) && ~isempty(index2)% if they are not empty
        opti.subject_to(kpe_var(index1)== kpe_var(index2));
        opti.subject_to(so_var(index1) == so_var(index2));
        opti.subject_to(sM_var(index1) == sM_var(index2));
        % we did not add additional constraints in lMo or lTs
        disp(['Constraint on muscles ' MuscleNames_total(index1) MuscleNames_total(index2)])
        end
    end
end
end
end
Sol = opti.solve();

%% Get and store solutions
clear Results

if strcmp(Misc.Mode_opti,'none')
    labelMRS='genericMRS';
elseif strcmp(Misc.Mode_opti,'passive')
    labelMRS='calibratedMRS';
    Results.bound.lMo =lMo_per_lim;
    Results.bound.lTs =lTs_per_lim;
    Results.bound.kpe =kpe_bound;
    Results.bound.so  =so_bound;
    Results.bound.sM  =sM_bound;
end

% Extract variables
lMo_calibrated = Sol.value(lMo_var);
lTs_calibrated = Sol.value(lTs_var);    
kpe_calibrated = Sol.value(kpe_var);
so_calibrated  = Sol.value(so_var);
sM_calibrated  = Sol.value(sM_var);

% Extract results
for trial=1:Misc.nTrials   
    % Store in a struct
    Results.MActivation(trial).(labelMRS)=Sol.value(TR(trial).a);
    Results.lMtildeopt(trial).(labelMRS) =Sol.value(TR(trial).lMtilde);
    Results.RMoment(trial).(labelMRS)    =Sol.value(TR(trial).aTk.*Misc.Topt);
    Results.MMoment(trial).(labelMRS)    =Sol.value(TR(trial).TsimVect);
    
    param_list=Sol.value(MuscProperties_params_all)';
end

import org.opensim.modeling.*
trial      = 4;
MuscleNames= DatStore(trial).MuscleNames;
nMuscles   = length(MuscleNames);
param_kT  =Misc.kT(:,Misc.idx_allMuscleList{trial}')';


Results.MuscleNames      =MuscleNames;
Results.kT.(labelMRS)    =param_kT;
Results.params.(labelMRS)=param_list;

save(fullfile(Misc.OutPath,[Misc.OutName '.mat']),'DatStore','Misc','Results')

% % create a model - not finished
% clear muscleParams
% for i=1:length(param_label)
%     muscleParams.([param_label{i}])=param_list(i,:);
% end
% muscleNames  = DatStore(trial).MuscleNames;
% modelPath    = char(Misc.model_path); 
% outPath      = Misc.OutPath;          
% newModelFile = strrep(Misc.model_name,'.osim','_cal.osim');
% CalParamsToOsim(muscleParams,muscleNames,modelPath,outPath,newModelFile); % previous ver: ParamsToOsim
% % test the compatibility of the model
% import org.opensim.modeling.*;
% model = Model(fullfile([outPath],newModelFile));
end
function [muscle_coupling_passive]= muscleCouplingPassive(model_conf,dof_name,side)
% Ankle joint
% {'soleus_r'} {'tib_post_r'} {'tib_ant_r'} are separated
% this group is plantar: 'per_brev' and 'per_long', and 'per_tert_r' is
% dorsi, thus per was not coupled
muscle_coupling.ankle_angle.model2392={
    'per_brev' 'per_long';
    'flex_dig' 'flex_hal';
    'ext_dig' 'ext_hal';
    'med_gas' 'lat_gas';
    };
muscle_coupling.ankle_angle.modelrajagopal={
    'perbrev' 'perlong'
    'fdl'     'fhl'
    'edl'     'ehl'
    'gasmed'  'gaslat'
};

% Knee joint
% {'sar_r'}    {'tfl_r'}    {'grac_r'}    {'rect_fem_r'} are separated
% this group is biarticular: 'semimem' 'semiten' 'bifemlh', and bifemsh is
% uniarticular, thus they were not coupled
muscle_coupling.knee_angle.model2392={
    'semimem' 'semiten';
    'semimem' 'bifemlh';
    'vas_med' 'vas_int';
    'vas_med' 'vas_lat';
    };
muscle_coupling.knee_angle.modelrajagopal={
    'semimem' 'semiten';
    'semimem' 'bflh';
    'vasmed'  'vasint';
    'vasmed'  'vaslat';
    };
% Hip joint
% this muscle group is joint: {'add_long_r'} {'add_brev_r'} is separated , and add_mag is too broad
% other muscles have too different muscle paths
muscle_coupling.hip_flexion.model2392={
    'glut_max1' 'glut_max2';
    'glut_max1' 'glut_max3';
    'glut_med1' 'glut_med2';
    'glut_med1' 'glut_med3';
    'glut_min1' 'glut_min2';
    'glut_min1' 'glut_min3';
    'add_mag1'  'add_mag2';
    'add_mag1'  'add_mag3';
    };
muscle_coupling.hip_flexion.modelrajagopal={
    'glmax1' 'glmax2';
    'glmax1' 'glmax3';
    'glmed1' 'glmed2';
    'glmed1' 'glmed3';
    'glmin1' 'glmin2';
    'glmin1' 'glmin3';
    'addmagDist'  'addmagIsch';
    'addmagDist'  'addmagMid';
    'addmagDist'  'addmagProx'; % one more than in 2392
    };

muscle_coupling.hip_adduction.model2392={};
muscle_coupling.hip_adduction.modelrajagopal={};
muscle_coupling.hip_rotation.model2392={};
muscle_coupling.hip_rotation.modelrajagopal={};

muscle_coupling_passive=append(muscle_coupling.(dof_name).(['model' model_conf{:}]),['_' side]); 
end