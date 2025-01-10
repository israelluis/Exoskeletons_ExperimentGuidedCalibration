function [Results,DatStore,Misc] = solveMuscleRedundancy_ExoCal(model_path,time,OutPath,Misc)
% -----------------------------------------------------------------------%
% INPUTS:
%           model_path: path to the .osim model
%           time: time window
%           OutPath: path to folder where results will be saved
%           Misc: structure with general input data (see manual for more details)
%
% OUTPUTS:
%           Results:    structure with outputs (states, controls, ...)
%           DatStore:   structure with data used for solving the optimal control problem
%           Misc:       structure with general input data (see manual for more details)
% -----------------------------------------------------------------------%


% update default settings
Misc = DefaultSettings_upd(Misc);

% number of motion trials
Misc.nTrials = length(Misc.IKfile);

%check if we have to adapt the start and end time so that it corresponds to
%time frames in the IK solution
[time] = Check_TimeIndices(Misc,time);
Misc.time=time;

% add the model name to the misc structure (mainly for post processing)
Misc.model_path = model_path;

%% Extract muscle information
% ----------------------------------------------------------------------- %
% Perform muscle analysis for the different selected trials
DatStore = struct;
for i = 1:Misc.nTrials
    % select the IK and ID file
    IK_path_trial = Misc.IKfile{i};
    ID_path_trial = Misc.IDfile{i};
    % Run muscle analysis
    % OPTION 0: Dont run and get the muscle analysis
    % OPTION 1: Run and get muscle analysis
    % OPTION 2: Dont run and get the muscle analysis from selected muscle analysis
    if Misc.GetAnalysis==0    
        OutPath_muscleAnalysis= OutPath;
        Misc.RunAnalysis      = 0;
    elseif Misc.GetAnalysis==1
        OutPath_muscleAnalysis= OutPath;
        Misc.RunAnalysis      = 1;
    elseif Misc.GetAnalysis==2
        Misc.RunAnalysis      = 0;
        dir_opt= [Misc.OutFolder '' Misc.OutFile(1:6) '_exoNN_simCAL_funACT_mMOO_\Results'];
        OutPath_muscleAnalysis= dir_opt;%obtain from SIMCAL
        if ~exist(OutPath,'dir'); mkdir(OutPath); end
    end
    MuscleAnalysisPath=fullfile(OutPath_muscleAnalysis,'MuscleAnalysis'); if ~exist(MuscleAnalysisPath,'dir'); mkdir(MuscleAnalysisPath); end
    if Misc.RunAnalysis
        disp('MuscleAnalysis Running .....');
        OpenSim_Muscle_Analysis(IK_path_trial,model_path,MuscleAnalysisPath,[time(i,1) time(i,end)],Misc.DofNames_Input)
        disp('MuscleAnalysis Finished');
    end
    Misc.MuscleAnalysisPath=MuscleAnalysisPath;
    
    % ----------------------------------------------------------------------- %
    % Extract muscle information -------------------------------------------- %
    % Get number of degrees of freedom (dofs), muscle-tendon lengths and moment
    % arms for the selected muscles.
    [~,Misc.trialName,~]=fileparts(IK_path_trial);
    if ~isfield(Misc,'MuscleNames_Input') || isempty(Misc.MuscleNames_Input)
        Misc=getMuscles4DOFS(Misc);
    end
    
    % get muscle geometry information
    [DatStore] = getMuscleInfo(IK_path_trial,ID_path_trial,Misc,DatStore,i);
    
    % display warnings in muscle selection
    Warnings_MuscleNames(DatStore,Misc,i);
    
    % get indexes of the muscles for which optimal fiber length, tendon stiffness are estimated
    [DatStore] = GetIndices_US(DatStore,Misc,i);
end

% set the tendon stiffness
if ~isfield(Misc,'kT') || isempty(Misc.kT)
    Misc.kT =ones(1,length(Misc.MuscleNames_Input)).*35;
end
if isfield(Misc,'Set_kT_ByName') && ~isempty(Misc.Set_kT_ByName)
    Misc= set_kT_ByName(Misc,DatStore);
end

% Shift tendon force-length curve as a function of the tendon stiffness
Misc.shift = getShift(Misc.kT);

% get the EMG information
[DatStore] = GetEMGInfo(Misc,DatStore);
[DatStore] = GetUSInfo(Misc,DatStore);

% get the number of muscles
NMuscles = length(DatStore(1).MuscleNames);

%% Static optimization
% ----------------------------------------------------------------------- %
% Solve the muscle redundancy problem using static optimization
% NOTE: We do not estimate any parameters here, but these results can serve as
% decent initial guess for the later dynamic optimization
% Extract the muscle-tendon properties
[Misc.params,Misc.lMo,Misc.lTs,Misc.FMo,Misc.alphao,Misc.kT,Misc.shift]=ReadMuscleParameters_upd(model_path,DatStore(1).MuscleNames,Misc);

[Misc]=ApplyCustomMadeForceChanges(Misc,DatStore(1).MuscleNames);

% Static optimization using IPOPT solver (used as an initial guess)
for trial = 1:Misc.nTrials
    DatStore    = SolveStaticOptimization_IPOPT_CasADi(DatStore,Misc,trial,0);
end

%% Input activation and contraction dynamics
% ----------------------------------------------------------------------- %
tau_act = 0.015;    Misc.tauAct = tau_act * ones(NMuscles, 1);       % activation time constant (activation dynamics)
tau_deact = 0.06;   Misc.tauDeact = tau_deact * ones(NMuscles,1);  % deactivation time constant (activation dynamics)
Misc.b = 0.1;       % tanh coefficient for smooth activation dynamics

%% Descretisation
Misc.Mesh_Frequency=100; %200
% mesh descretisation
for trial = 1:Misc.nTrials
    t0 = DatStore(trial).time(1); tf = DatStore(trial).time(end);
    Mesh(trial).N = round((tf-t0)*Misc.Mesh_Frequency);
    Mesh(trial).step = (tf-t0)/Mesh(trial).N;
    Mesh(trial).t = t0:Mesh(trial).step:tf;
end

%% Evaluate splines at Mesh Points
% ----------------------------------------------------------------------- %
% Get IK, ID, muscle analysis and static opt information at mesh points

for trial = 1:Misc.nTrials
    % Discretization
    N = Mesh(trial).N;
    time_opt = Mesh(trial).t;
    % Spline approximation of muscle-tendon length (LMT), moment arms (MA) and inverse dynamic torques (ID)
    for dof = 1:DatStore(trial).nDOF
        for m = 1:NMuscles
            DatStore(trial).JointMASpline(dof).Muscle(m) = spline(DatStore(trial).time,squeeze(DatStore(trial).dM(:,dof,m)));
        end
        DatStore(trial).JointIDSpline(dof) = spline(DatStore(trial).time,DatStore(trial).T_exp(:,dof));
    end
    
    for m = 1:NMuscles
        DatStore(trial).LMTSpline(m) = spline(DatStore(trial).time,DatStore(trial).LMT(:,m));
    end
    
    % Evaluate LMT, VMT, MA and ID at optimization mesh
    DatStore(trial).LMTinterp = zeros(length(time_opt),NMuscles); % Muscle-tendon length
    for m = 1:NMuscles
        [DatStore(trial).LMTinterp(:,m),~,~] = SplineEval_ppuval(DatStore(trial).LMTSpline(m),time_opt,1);
    end
    DatStore(trial).MAinterp = zeros(length(time_opt),DatStore(trial).nDOF*NMuscles); % Moment arm
    DatStore(trial).IDinterp = zeros(length(time_opt),DatStore(trial).nDOF); % Inverse dynamic torque
    for dof = 1:DatStore(trial).nDOF
        for m = 1:NMuscles
            index_sel=(dof-1)*NMuscles+m;
            DatStore(trial).MAinterp(:,index_sel) = ppval(DatStore(trial).JointMASpline(dof).Muscle(m),time_opt);
        end
        DatStore(trial).IDinterp(:,dof) = ppval(DatStore(trial).JointIDSpline(dof),time_opt);
    end
    
    % Interpolate results of static optimization
    DatStore(trial).SoActInterp = interp1(DatStore(trial).time,DatStore(trial).SoAct,time_opt');
    DatStore(trial).SoRActInterp = interp1(DatStore(trial).time,DatStore(trial).SoRAct,time_opt');
    DatStore(trial).SoForceInterp = interp1(DatStore(trial).time,DatStore(trial).SoForce.*DatStore(trial).cos_alpha./Misc.FMo,time_opt);
    [~,DatStore(trial).lMtildeInterp ] = FiberLength_Ftilde(DatStore(trial).SoForceInterp,Misc.params,DatStore(trial).LMTinterp,Misc.kT,Misc.shift);
    DatStore(trial).vMtildeinterp = zeros(size(DatStore(trial).lMtildeInterp));
    for m = 1:NMuscles
        DatStore(trial).lMtildeSpline = spline(time_opt,DatStore(trial).lMtildeInterp(:,m));
        [~,DatStore(trial).vMtildeinterp_norm,~] = SplineEval_ppuval(DatStore(trial).lMtildeSpline,time_opt,1);
        DatStore(trial).vMtildeinterp(:,m) = DatStore(trial).vMtildeinterp_norm;
    end
    
    % interpolate the joint angles
    DatStore(trial).IKinterp = interp1(DatStore(trial).time,DatStore(trial).q_exp,time_opt');
end

%% Fiber length setup
% Obtain index and interpolate for muscle fiber calibration
if Misc.muscleFiberCal==1
    [US_data,nUSdata,ind_US,ind_US_AT,ind_USnone,USDigitalizedInterp]=fiberCalibrationSetup(DatStore,Misc,time_opt); % only works for one trial
end
%% setup options for the solver
% Create an NLP solver
% output.setup.lM_projecteddata = lM_projecteddata;
output.setup.nlp.solver = 'ipopt';
output.setup.nlp.ipoptoptions.linear_solver = 'mumps';
% Set derivativelevel to 'first' for approximating the Hessian
output.setup.derivatives.derivativelevel = 'second';
output.setup.nlp.ipoptoptions.tolerance = 1e-6;
output.setup.nlp.ipoptoptions.maxiterations = 1000;
if strcmp(output.setup.derivatives.derivativelevel, 'first')
    optionssol.ipopt.hessian_approximation = 'limited-memory';
end
% By default, the barrier parameter update strategy is monotone.
% https://www.coin-or.org/Ipopt/documentation/node46.html#SECTION000116020000000000000
% Uncomment the following line to use an adaptive strategy
% optionssol.ipopt.mu_strategy = 'adaptive'; 
optionssol.ipopt.nlp_scaling_method = 'gradient-based';
optionssol.ipopt.linear_solver = output.setup.nlp.ipoptoptions.linear_solver;
optionssol.ipopt.tol = output.setup.nlp.ipoptoptions.tolerance;
optionssol.ipopt.max_iter = output.setup.nlp.ipoptoptions.maxiterations;

%% Dynamic Optimization - Default parameters
% ----------------------------------------------------------------------- %
% Solve muscle redundancy problem with default parameters
% Problem bounds
e_min = 0; e_max = 1;                   % bounds on muscle excitation
a_min = 0; a_max = 1;                   % bounds on muscle activation
vMtilde_min = -10; vMtilde_max = 10;    % bounds on normalized muscle fiber velocity
lMtilde_min = 0.1; lMtilde_max = 1.7;   % bounds on normalized muscle fiber length


% CasADi setup
import casadi.*
opti    = casadi.Opti();    % create opti structure

% get total number of mesh points
nTrials = Misc.nTrials;
N_tot = sum([Mesh().N]);

Misc.Topt=100; % 100 (150 byDefault)

% set intial guess based on static opt data
SoActGuess = zeros(NMuscles,N_tot);
SoExcGuess = zeros(NMuscles,N_tot-nTrials);
lMtildeGuess = zeros(NMuscles,N_tot);
vMtildeGuess = zeros(NMuscles,N_tot-nTrials);
SoRActGuess = zeros(DatStore(1).nDOF,N_tot-nTrials);
ctx = 1;  ctu= 1;

% if Misc.Exo_Enable==0; IG_source='SO'; else; IG_source='noExo'; end

IG_source='SO'; %noExo SO
if strcmp(IG_source,'noExo')
    dir_opt= [Misc.OutFolder '\' Misc.OutFile(1:6) '_exoNN_simCAL_funACT_ORun\Results']; %exoNN_simCAL_funACT
    store_name='simulation_GenCalValResults';
    sim=load(fullfile(dir_opt,store_name));
for trial = 1:nTrials
    ctx_e = ctx+Mesh(trial).N;      % counter for states
    ctu_e = ctu+Mesh(trial).N-1;    % counter for controls
    SoActGuess(:,ctx:ctx_e) = sim.Results.MActivation.genericMRS;
    SoExcGuess(:,ctu:ctu_e) = sim.Results.MActivation.genericMRS(:,1:end-1);
    lMtildeGuess(:,ctx:ctx_e) = sim.Results.lMtildeopt.genericMRS;
    vMtildeGuess(:,ctu:ctu_e) = sim.Results.vMtilde.genericMRS;
    SoRActGuess(:,ctu:ctu_e) = sim.Results.RActivation.genericMRS; %RA torque, not activation
    ctx = ctx_e+1;
    ctu = ctu_e+1;
end    
elseif strcmp(IG_source,'SO')
for trial = 1:nTrials
    ctx_e = ctx+Mesh(trial).N;      % counter for states
    ctu_e = ctu+Mesh(trial).N-1;    % counter for controls
    SoActGuess(:,ctx:ctx_e) = DatStore(trial).SoActInterp';
    SoExcGuess(:,ctu:ctu_e) = DatStore(trial).SoActInterp(1:end-1,:)';
    lMtildeGuess(:,ctx:ctx_e) = DatStore(trial).lMtildeInterp';
    vMtildeGuess(:,ctu:ctu_e) = DatStore(trial).vMtildeinterp(1:end-1,:)';
    SoRActGuess(:,ctu:ctu_e) = DatStore(trial).SoRActInterp(1:end-1,:)';
    ctx = ctx_e+1;
    ctu = ctu_e+1;
end
end

% States
%   - muscle activation
a = opti.variable(NMuscles,N_tot+nTrials);      % Variable at mesh points
opti.subject_to(a_min < a < a_max);             % Bounds
opti.set_initial(a,SoActGuess);                 % Initial guess (static optimization)
%   - Muscle fiber lengths
lMtilde = opti.variable(NMuscles,N_tot+nTrials);
opti.subject_to(lMtilde_min < lMtilde < lMtilde_max);
opti.set_initial(lMtilde,lMtildeGuess);
%   - Controls
e = opti.variable(NMuscles,N_tot);
opti.subject_to(e_min < e < e_max);
opti.set_initial(e, SoExcGuess);
%   - Reserve actuators
aT = opti.variable(DatStore(trial).nDOF,N_tot);
opti.subject_to(-1 < aT <1);
% opti.set_initial(aT, SoRActGuess/Misc.Topt); %Not added. Longer to converge
%   - Time derivative of muscle-tendon forces (states)
vMtilde = opti.variable(NMuscles,N_tot);
opti.subject_to(vMtilde_min < vMtilde < vMtilde_max);
opti.set_initial(vMtilde,vMtildeGuess);
%   - Projected muscle fiber length - Auxilary variable to avoid muscle buckling & square root expression in muscle dynamics
lM_projected = opti.variable(NMuscles,N_tot + nTrials);
opti.subject_to(1e-4 < lM_projected(:)); % We impose that projected muscle fiber length has strict positive length

%% Declare and setup muscle fiber optimization
J=0;
if Misc.muscleFiberCal==1
    [opti,J,kT,shift,lMo,lTs] = formulateMuscleFiberCalibration(opti,J,Misc,lMtilde,nUSdata,ind_US,ind_USnone,USDigitalizedInterp,N,NMuscles,DatStore.MuscleNames);
else
    [opti,J,kT,shift,lMo,lTs] = formulateMuscleFiberCalibration(opti,J,Misc,lMtilde,[],[],[],[],[],[],[]);
end

%% Initial guess for this variable is retrieved from lMtilde guess
% and geometric relationship between pennation angle, muscle length
% and width
lMo_default = Misc.params(2,:)'; 
alphao = Misc.params(4,:)';
lMGuess = lMtildeGuess.*lMo_default;
w = lMo_default.*sin(alphao);
lM_projectedGuess = sqrt((lMGuess.^2 - w.^2));
opti.set_initial(lM_projected,lM_projectedGuess);

% constraint on projected fiber length
w = lMo.*sin(alphao);
lM = lMtilde.*lMo;
opti.subject_to(lM.^2 - w.^2 == lM_projected.^2);

% output optimization variables
MuscProperties_params =[Misc.params(1,:)' lMo lTs Misc.params(4:8,:)']';
MuscProperties_kT     =kT';
MuscProperties_shift  =shift';
%% Exoskeleton modelling
if Misc.Exo_Enable
    for trial = 1:Misc.nTrials
        NExo=length(Misc.Exo_Dof);
                
        % get kinematics and time
        iSel = find(ismember(Misc.DofNames_Input,Misc.Exo_Dof));
        qSel = DatStore(trial).IKinterp(1:N_tot,iSel);
        time_mesh=Mesh(trial).t(1:N_tot);
        t_toe_off=Misc.toe_off_timing;
        
        % initialize variables
        if strcmp(Misc.Exo_Type,'active') && strcmp(Misc.Exo_Mode,'optimal') 
            idealT = opti.variable(NExo,N_tot);
        elseif strcmp(Misc.Exo_Type,'active') && strcmp(Misc.Exo_Mode,'manual') 
            idealT = NaN(NExo,N_tot);
        elseif strcmp(Misc.Exo_Type,'quasi') && strcmp(Misc.Exo_Mode,'optimal') 
            Exo_kvalue= opti.variable(trial,NExo);  %stiffness
            timeExo_i = opti.variable(trial,NExo);  %engage clutching time
            timeExo_f = opti.variable(trial,NExo);  %disengage clutching time   
            
            Exo_kvalue_lim=NaN(NExo,2);
            timeExo_i_lim =NaN(NExo,2);
            timeExo_f_lim =NaN(NExo,2);

            Exo_kvalue_IG =NaN(NExo,1);
            timeExo_i_IG  =NaN(NExo,1);
            timeExo_f_IG  =NaN(NExo,2);                                                                                        % it should be one! 
        elseif strcmp(Misc.Exo_Type,'quasi') && (strcmp(Misc.Exo_Mode,'manual') || strcmp(Misc.Exo_Mode,'informed')) 
            Exo_kvalue= NaN(trial,NExo);            %stiffness
            timeExo_i = NaN(trial,NExo);            %engage clutching time
            timeExo_f = NaN(trial,NExo);            %disengage clutching time         
        end
        
        % loop for each exoskeleton
        for ExoSel=1:length(Misc.Exo_Dof)
            Exo_Dof_sel=Misc.Exo_Dof(ExoSel);    
            group_sign =Misc.Exo_group(ExoSel);
            
            % configuration for exoskeleton
            if strcmp(Misc.Exo_Type,'active') && strcmp(Misc.Exo_Mode,'manual') 
                torque=Misc.torqueProfile{ExoSel};
                %interpolate to mesh points
                idx            = 1:length(torque);                                 
                idxq           = linspace(min(idx), max(idx), N_tot); 
                Vi             = interp1(idx, torque, idxq, 'pchip');
            	idealT(ExoSel,:)=Vi;
            elseif  strcmp(Misc.Exo_Type,'active') && strcmp(Misc.Exo_Mode,'optimal')  
                range_exo=[0 250]*group_sign;
                range_exo = sort(range_exo);
                opti.subject_to(min(range_exo) < idealT(ExoSel,:) < max(range_exo));
            elseif strcmp(Misc.Exo_Type,'quasi') && strcmp(Misc.Exo_Mode,'informed')    
                % INPUT : k, time_i and time_f optimized
                % OUTPUT: qVal @ time_i and time_f
                Exo_kvalue(trial,ExoSel) = Misc.Exo_kvalue(ExoSel);
                timeExo_i(trial,ExoSel)  = Misc.Exo_time_i(ExoSel);
                timeExo_f(trial,ExoSel)  = Misc.Exo_time_f(ExoSel);
                
                qSel_val = qSel(:,ExoSel);  
                pp       = spline(time_mesh,qSel_val);
                qOffExo_i= ppval(pp,timeExo_i(trial,ExoSel));
                qOffExo_f= ppval(pp,timeExo_f(trial,ExoSel));
                
            elseif strcmp(Misc.Exo_Type,'quasi') && strcmp(Misc.Exo_Mode,'manual')
                % INPUT : k and time_i manual
                % OUTPUT: time_f, and qVal @ time_i and time_f
                Exo_kvalue(trial,ExoSel) = Misc.Exo_kvalue(ExoSel);
                timeExo_i                = Misc.Exo_time(ExoSel);
                qSel_val                 = qSel(:,ExoSel);       
%               we could use spline but then it has to be inverted with slmsolve
%               to find the x value at which the angle is zero. Instead I approximate
                frame_i   =find(time_mesh>=timeExo_i,1,'first');
                qOffExo_i =qSel_val(frame_i);
                timeExo_i(trial,ExoSel) =time_mesh(frame_i);
                
                frame_f   =find((-group_sign.*qSel_val(frame_i:end))<-group_sign.*qOffExo_i,1,'first');
                if isempty(frame_f) % if there is no value
                    frame_f= frame_i; 
                else                % if there is value
                    frame_f= frame_i+frame_f-2;
                end
                qOffExo_f=qSel_val(frame_f);
                timeExo_f(trial,ExoSel) =time_mesh(frame_f);  
                
            elseif strcmp(Misc.Exo_Type,'quasi') && strcmp(Misc.Exo_Mode,'optimal') 
                if strcmp(Exo_Dof_sel{1}(1:end-2),'ankle_angle') % strategy for ankle plantarflexors
%                     plot(time_mesh,qSel); hold on; plot(time_mesh(ind_min),qSel(ind_min),'or'); plot(time_mesh(ind_max),qSel(ind_max),'or'); xline(t_toe_off)
                    % get angles
                    [~,ind_max]=max(qSel(:,ExoSel));
                    [~,ind_min]=min(qSel(1:ind_max,ExoSel));
                    t_min_ank =time_mesh(ind_min);
                    t_max_ank =time_mesh(ind_max);
                    % kvalue
                    Exo_kvalue_lim(ExoSel,:)=[0 400];
                    Exo_kvalue_IG(ExoSel)   =200;
                    % engage
                    timeExo_i_lim(ExoSel,:) =[t_min_ank t_max_ank];
                    timeExo_i_IG(ExoSel)    =t_min_ank;
                    % disengage
                    timeExo_f_lim(ExoSel,:) =[t_max_ank t_toe_off];
                    timeExo_f_IG(ExoSel)    =mean([t_max_ank t_toe_off]);

                elseif strcmp(Exo_Dof_sel{1}(1:end-2),'hip_flexion') % strategy for hip flexors   
%                     plot(time_mesh,qSel); hold on; xline(mean([time_mesh(1) t_toe_off])); xline(t_toe_off)
                    [~,ind_toe_off]= min(abs(time_mesh-t_toe_off));
                    [~,ind_max1]   = max(qSel(1:ind_toe_off  ,ExoSel));
                    [~,ind_max2]   = max(qSel(ind_toe_off:end,ExoSel));
                    
                    [~,ind_minPeak]  =min(qSel(1:ind_toe_off  ,ExoSel));
                    
                    % kvalue
                    Exo_kvalue_lim(ExoSel,:)=[0 400];
                    Exo_kvalue_IG(ExoSel)   =100;
                    % engage
                    timeExo_i_lim(ExoSel,:) =[time_mesh(ind_max1) time_mesh(ind_minPeak)];
                    timeExo_i_IG(ExoSel)    =mean([time_mesh(ind_max1) time_mesh(ind_minPeak)]);
                    % disengage
                    timeExo_f_lim(ExoSel,:) =[time_mesh(ind_minPeak) time_mesh(ind_toe_off+ind_max2-1)];
                    timeExo_f_IG(ExoSel)    = t_toe_off; % you can try algo: mean([t_toe_off time_mesh(end)])           
                    
                elseif strcmp(Exo_Dof_sel{1}(1:end-2),'knee_angle') % strategy for knee extensors
                    % get angles
                    ind_maxes=find(islocalmax(qSel(:,ExoSel)));
                    ind_mines=find(islocalmin(qSel(:,ExoSel)));
                    t_maxPeak1=time_mesh(ind_maxes(1));
                    t_minAng1 =time_mesh(ind_mines(1));
                    % kvalue
                    Exo_kvalue_lim(ExoSel,:)=[0 400];
                    Exo_kvalue_IG(ExoSel)   =200;
                    % engage
                    timeExo_i_lim(ExoSel,:) =[time_mesh(1) t_maxPeak1];
                    timeExo_i_IG(ExoSel)    = time_mesh(1);
                    % disengage
                    timeExo_f_lim(ExoSel,:) =[t_maxPeak1 t_toe_off];
                    timeExo_f_IG(ExoSel)    =t_minAng1;
                elseif strcmp(Exo_Dof_sel{1}(1:end-2),'hip_adduction') % strategy for hip abductors
                    % get angles
                    [~,ind_max]=max(qSel(:,ExoSel));
                    t_max_hip =time_mesh(ind_max);
                    % kvalue
                    Exo_kvalue_lim(ExoSel,:)=[0 500];
                    Exo_kvalue_IG(ExoSel)   =200;
                    % engage
                    timeExo_i_lim(ExoSel,:) =[time_mesh(1) t_max_hip];
                    timeExo_i_IG(ExoSel)    = time_mesh(1);
                    % disengage
                    timeExo_f_lim(ExoSel,:) =[t_max_hip t_toe_off];
                    timeExo_f_IG(ExoSel)    =t_toe_off;    
                end
                
               % apply limits
                opti.subject_to(Exo_kvalue_lim(ExoSel,1) < Exo_kvalue(trial,ExoSel) < Exo_kvalue_lim(ExoSel,end));
               
                opti.subject_to(timeExo_i_lim(ExoSel,1) < timeExo_i(trial,ExoSel) < timeExo_i_lim(ExoSel,end));
                opti.set_initial(timeExo_i(trial,ExoSel), timeExo_i_IG(ExoSel));

                opti.subject_to(timeExo_f_lim(ExoSel,1) < timeExo_f(trial,ExoSel) < timeExo_f_lim(ExoSel,end));
                opti.set_initial(timeExo_f(trial,ExoSel), timeExo_f_IG(ExoSel));
                    
                opti.subject_to(timeExo_f(trial,ExoSel) - timeExo_i(trial,ExoSel) > 0);
                
                % bspline to relate angle and time
                qSel_val  = qSel(:,ExoSel);
                pp_casadi = casadi.interpolant('pp', 'bspline', {time_mesh}, qSel_val);
                qOffExo_i = pp_casadi(timeExo_i(trial,ExoSel));
                qOffExo_f = pp_casadi(timeExo_f(trial,ExoSel));

                % engage and disengage at the same angle 
                opti.subject_to((qOffExo_i - qOffExo_f).^2 < 1e-1); %1e-2
            end
  
            if strcmp(Misc.Exo_Type,'active')
                Texo(trial,ExoSel).T          = idealT(ExoSel,:);
                Texo(trial,ExoSel).Ind        = iSel(ExoSel); % index of applied torque             
            elseif strcmp(Misc.Exo_Type,'quasi')
                % compute torque based on angle value clutching 
                b_exo= 100;
                BoolTime_high=tanh(b_exo*(time_mesh-timeExo_i(trial,ExoSel)))*0.5+0.5;
                BoolTime_low =tanh(b_exo*-(time_mesh-timeExo_f(trial,ExoSel)))*0.5+0.5;
                BoolTimeActive=BoolTime_high.*BoolTime_low;
                BoolqL = 1; 
%                 BoolqL =(tanh(1000*(-group_sign.*(qSel_val-qOffExo_i)))*0.5+0.5); % this is need for passive-only exoskeleton
                BoolActive=BoolqL.*BoolTimeActive';
    
                angle=-group_sign.*(qSel_val-qOffExo_i)*pi/180;
                angle_pos=-group_sign.*(tanh(b_exo*(angle))*0.5+0.5).*angle;
                Texo(trial,ExoSel).T          = -((Exo_kvalue(trial,ExoSel).*angle_pos).*BoolActive)';

                Texo(trial,ExoSel).Ind        = iSel(ExoSel); 
                Texo(trial,ExoSel).BoolActive = BoolActive;
                Texo(trial,ExoSel).qOffExo_i  = qOffExo_i;
                Texo(trial,ExoSel).qOffExo_f  = qOffExo_f;
                Texo(trial,ExoSel).qSel       = qSel_val;
            end
        end
    end
end
%% Muscle Factors
MuscProperties_factor=Misc.forGen.actFM.out.percentages;

%% Implemetation of controls, states and states derivatives
N_acc = 0; % Index that keeps track of trials that are accumulated
% Loop over trials --> one simulation for each trial
for trial = 1:Misc.nTrials
    % Time bounds
    
    t0 = DatStore(trial).time(1); tf = DatStore(trial).time(end);
    % Discretization
    N = Mesh(trial).N;
    h = Mesh(trial).step;
    
    % Loop over mesh points formulating NLP
    for k=1:N
        % Variables within current mesh interval
        ak = a(:,(N_acc+trial-1) + k); lMtildek = lMtilde(:,(N_acc+trial-1) + k);
        vMtildek = vMtilde(:,N_acc + k); aTk = aT(:,N_acc + k); ek = e(:,N_acc + k);
        lM_projectedk = lM_projected(:,(N_acc+trial-1) + k);
        
        % Euler integration  Uk = (X_(k+1) - X_k)/*dt
        Xk = [ak; lMtildek];
        Zk = [a(:,(N_acc+trial-1) + k + 1);lMtilde(:,(N_acc+trial-1) + k + 1)];
        Uk = [ActivationDynamics(ek,ak,Misc.tauAct,Misc.tauDeact,Misc.b); vMtildek];
        opti.subject_to(eulerIntegrator(Xk,Zk,Uk,h) == 0);
        
        % Get muscle-tendon forces and derive Hill-equilibrium
        [Hilldiffk,FTk] = ForceEquilibrium_lMtildeState_optPassive(ak,lMtildek,vMtildek,lM_projectedk,...
            DatStore(trial).LMTinterp(k,:)',MuscProperties_params',MuscProperties_factor,MuscProperties_kT',MuscProperties_shift'); 
        
        % Add path constraints
        % Moment constraints
        for dof = 1:DatStore(trial).nDOF
            T_exp = DatStore(trial).IDinterp(k,dof);
            index_sel = (dof-1)*(NMuscles)+1:(dof*NMuscles); % moment is a vector with the different dofs "below" each other
            T_sim = DatStore(trial).MAinterp(k,index_sel)*FTk + Misc.Topt*aTk(dof);
            
            % subtract exoskeleton moment from ID torque
            % (TID = Tmuscles + Texo)
            if Misc.Exo_Enable
                for ExoSel= 1:length(Misc.Exo_Dof)
                    if dof == Texo(trial,ExoSel).Ind
                        T_exp = T_exp - Texo(trial,ExoSel).T(k);
                    end
                end
            end
            opti.subject_to(T_exp - T_sim == 0);
        end
        % Hill-equilibrium constraint
        opti.subject_to(Hilldiffk == 0);
    end
    N_acc = N_acc + N;
% Cost function
J = J + ... 
    Misc.wAct*0.5*(sumsqr(e)/N/NMuscles + sumsqr(a)/N/NMuscles) + ...
    Misc.wTres*sumsqr(aT)/N/DatStore(trial).nDOF + ...
    Misc.wVm*sumsqr(vMtilde)/N/NMuscles; % this is faster than sumsqr(e) and sumsqr(a) has excitations fast spikes
end

Misc.wExo=0; %0.0001 % in case regularization is needed
if Misc.Exo_Enable && strcmp(Misc.Exo_Type,'quasi')
    exoR_acc=0;
    for trial = 1:Misc.nTrials
        time_mesh=Mesh(trial).t;
        for dof = 1:DatStore(trial).nDOF
            for ExoSel= 1:length(Misc.Exo_Dof)
                if dof == Texo(trial,ExoSel).Ind
                    gcPerExo=(timeExo_i(trial,ExoSel))./(time_mesh(end)-time_mesh(1));
                    exoR_acc=exoR_acc+gcPerExo;
                end
            end
        end
    end
    J = J + Misc.wExo*exoR_acc;
end

opti.minimize(J); % Define cost function in opti

% Create an NLP solver
opti.solver(output.setup.nlp.solver,optionssol);

% Solve
diary(fullfile(OutPath,[Misc.OutName 'GenericMRS.txt']));
tic
sol = opti.solve();
dt = toc;
disp(['Computation time solving OCP: ' num2str(dt) ' s'])
diary off

% Extract results
% Variables at mesh points
% Muscle activations and muscle-tendon forces
a_opt = sol.value(a);
lMtilde_opt = sol.value(lMtilde);
% Muscle excitations
e_opt = sol.value(e);
% Reserve actuators
aT_opt = sol.value(aT);
% Time derivatives of muscle-tendon forces
vMtilde_opt = sol.value(vMtilde);
% Optimal lM_projectedilary variable
lM_projected_opt = sol.value(lM_projected);

% get parameters
MuscProperties_params_opt=sol.value(MuscProperties_params);
MuscProperties_kT_opt    =sol.value(MuscProperties_kT);
MuscProperties_shift_opt =sol.value(MuscProperties_shift);

% store parameters
Results.params.calibratedMRS=MuscProperties_params_opt;
Results.kT.calibratedMRS    =MuscProperties_kT_opt;
Results.params.genericMRS   =Misc.params;
Results.kT.genericMRS       =Misc.kT;

% compute unNormalized tendon values
if Misc.muscleFiberCal == 1
kT_sel =MuscProperties_kT_opt(ind_US_AT);
FMo_sel=MuscProperties_params_opt(1,ind_US_AT);
lTs_sel=MuscProperties_params_opt(3,ind_US_AT);
[kT_abs_computed] = unNormalizedTendon(kT_sel,FMo_sel,lTs_sel);
Results.kT_ATunNormalized.calibratedMRS=sol.value(kT_abs_computed);

kT_sel =Misc.kT(ind_US_AT);
FMo_sel=Misc.params(1,ind_US_AT);
lTs_sel=Misc.params(3,ind_US_AT);
[kT_abs_computed] = unNormalizedTendon(kT_sel,FMo_sel,lTs_sel);
Results.kT_ATunNormalized.genericMRS   =sol.value(kT_abs_computed);
% are they differents? quick to check -> compare_param=MuscProperties_params_opt-Misc.params;
end

% Append results to output structures of exoskeleton
if Misc.Exo_Enable
    for trial = 1:Misc.nTrials
        for ExoSel= 1:length(Misc.Exo_Dof)
            if strcmp(Misc.Exo_Type,'active') && (strcmp(Misc.Exo_Mode,'optimal') || strcmp(Misc.Exo_Mode,'manual'))
                % no parameters related to timing, angle or stiffness
                % torque
                Results.Texo(trial,ExoSel).T          = sol.value(Texo(trial,ExoSel).T);   
            elseif strcmp(Misc.Exo_Type,'quasi') && strcmp(Misc.Exo_Mode,'optimal')
                % parameters
                Results.Exo(trial,ExoSel).t_i =sol.value(timeExo_i(trial,ExoSel));
                Results.Exo(trial,ExoSel).t_f =sol.value(timeExo_f(trial,ExoSel));        
                Results.Exo(trial,ExoSel).q_i =sol.value(Texo(trial,ExoSel).qOffExo_i);
                Results.Exo(trial,ExoSel).q_f =sol.value(Texo(trial,ExoSel).qOffExo_f);
                Results.Exo(trial,ExoSel).qSel=Texo(trial,ExoSel).qSel ;    
                Results.Exo(trial,ExoSel).k   =sol.value(Exo_kvalue(trial,ExoSel));    
                % torque and boolean
                Results.Texo(trial,ExoSel).T          = sol.value(Texo(trial,ExoSel).T);
                Results.Texo(trial,ExoSel).BoolActive = sol.value(Texo(trial,ExoSel).BoolActive);   
            elseif strcmp(Misc.Exo_Type,'quasi') && (strcmp(Misc.Exo_Mode,'manual') || strcmp(Misc.Exo_Mode,'informed')) 
                % parameters
                Results.Exo(trial,ExoSel).t_i =timeExo_i(trial,ExoSel);
                Results.Exo(trial,ExoSel).t_f =timeExo_f(trial,ExoSel);        
                Results.Exo(trial,ExoSel).q_i =Texo(trial,ExoSel).qOffExo_i;
                Results.Exo(trial,ExoSel).q_f =Texo(trial,ExoSel).qOffExo_f;
                Results.Exo(trial,ExoSel).qSel=Texo(trial,ExoSel).qSel ;    
                Results.Exo(trial,ExoSel).k   =Exo_kvalue(trial,ExoSel);      
                % torque and boolean
                Results.Texo(trial,ExoSel).T          = Texo(trial,ExoSel).T;
                Results.Texo(trial,ExoSel).BoolActive = Texo(trial,ExoSel).BoolActive;     
            end
        end
    end
end

% Append results to output structures
Ntot = 0;
for trial = 1:nTrials
    t0 = DatStore(trial).time(1); tf = DatStore(trial).time(end);
    N = round((tf-t0)*Misc.Mesh_Frequency);
    % Time grid
    tgrid = linspace(t0,tf,N+1)';
    % Save results
    Results.Time(trial).genericMRS = tgrid;
    Results.MActivation(trial).genericMRS = a_opt(:,(Ntot + trial - 1) + 1:(Ntot + trial - 1) + N + 1);
    Results.lMtildeopt(trial).genericMRS = lMtilde_opt(:,(Ntot + trial - 1) + 1:(Ntot + trial - 1) + N + 1);
    Results.lM(trial).genericMRS = lMtilde_opt(:,(Ntot + trial - 1) + 1:(Ntot + trial - 1) + N + 1).*repmat(Misc.lMo',1,length(tgrid));
    Results.vMtilde(trial).genericMRS = vMtilde_opt(:,Ntot + 1:Ntot + N);
    Results.MExcitation(trial).genericMRS = e_opt(:,Ntot + 1:Ntot + N);
    Results.RActivation(trial).genericMRS = aT_opt(:,Ntot + 1:Ntot + N)*Misc.Topt;
    Results.MuscleNames = DatStore.MuscleNames;
    Results.OptInfo = output;
    % Tendon force
    Results.lMTinterp(trial).genericMRS = DatStore(trial).LMTinterp';
    [TForcetilde_,TForce_] = TendonForce_lMtilde(Results.lMtildeopt(trial).genericMRS',MuscProperties_params_opt,Results.lMTinterp(trial).genericMRS',MuscProperties_kT_opt,MuscProperties_shift_opt);
    Results.TForcetilde(trial).genericMRS = TForcetilde_';
    Results.TForce(trial).genericMRS = TForce_';
    % get information F/l and F/v properties - updated with passive forces
    [Fpe_,FMltilde_,FMvtilde_] = getForceLengthVelocityProperties_setPassiveParam(Results.lMtildeopt(trial).genericMRS',Results.vMtilde(trial).genericMRS',MuscProperties_params_opt(5,:),...
                                                                                  MuscProperties_params_opt(6,:),MuscProperties_params_opt(7,:),MuscProperties_params_opt(8,:));
    FMo = ones(N+1,1)*Misc.params(1,:);
    Results.Fpe(trial).genericMRS = (Fpe_.*FMo)';
    Results.FMltilde(trial).genericMRS = FMltilde_';
    Results.FMvtilde(trial).genericMRS = FMvtilde_';
    Ntot = Ntot + N;
end

%% Store Results

% store the Misc structure as well in the results
Results.Misc = Misc;

% add selected muscle names to the output structure
Results.MuscleNames = DatStore.MuscleNames;

%% save the results
% plot states and variables from parameter estimation simulation
save(fullfile(OutPath,[Misc.OutName 'Results.mat']),'Results','DatStore','Misc');


end


function [Misc] = ApplyCustomMadeForceChanges(Misc,names)

names_upd =Misc.forGen.actFM.names;
% update if all muscle include
if strcmp(names_upd,'all');names_upd= cellfun(@(x) x(1:end-2), names, 'un', 0); end
nMuscles_upd  =length(names_upd);
nMuscles_all  =length(names);

ind=zeros(nMuscles_upd,1);
for j=1:nMuscles_upd
    ind(j)=find(strcmp(names,[names_upd{j} '_' Misc.side_sel])); %
end

Misc.forGen.actFM.out.percentages=ones(1,nMuscles_all);
Misc.forGen.actFM.out.percentages(ind)=Misc.forGen.actFM.percentages;
end