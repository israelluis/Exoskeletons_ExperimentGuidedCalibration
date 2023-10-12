function [opti,J,kT,shift,lMo,lTs] = formulateMuscleFiberCalibration(opti,J,Misc,lMtilde,nUSdata,ind_US,ind_USnone,USDigitalizedInterp,N,NMuscles,MuscleNames)
import casadi.*

% read variables
kT_list   = Misc.kT';
shift_list= Misc.shift';
% FMo_list  = Misc.params(1,:)';
lMo_list  = Misc.params(2,:)';
lTs_list  = Misc.params(3,:)';

% set optimization variables
kT        =MX(40,1);
shift     =MX(40,1);
lMo       =MX(40,1);
lTs       =MX(40,1);
if Misc.muscleFiberCal == 0
    kT    =kT_list;
    shift =shift_list;    
    lMo   =lMo_list;
    lTs   =lTs_list;    
elseif Misc.muscleFiberCal == 1
    Misc.wUS =  1.0; %  Misc.wUS =0.1:minEffort is stronger (82 sec)   =1:good (82 sec)   =10:marginally-better (102 sec)
    kT_opt  = opti.variable(nUSdata,1);
    lMo_opt = opti.variable(nUSdata,1);
    lTs_opt = opti.variable(nUSdata,1);
    
    kT(ind_US)     =kT_opt;
    kT(ind_USnone) =kT_list(ind_USnone);
    shift          =getShift(kT); 
    
    lMo(ind_US)    =lMo_opt;
    lMo(ind_USnone)=lMo_list(ind_USnone);
    lTs(ind_US)    =lTs_opt;
    lTs(ind_USnone)=lTs_list(ind_USnone);
    
    % ward et al.
    MuscleNames_sel = MuscleNames(ind_US);
    lMo_ward        = getOptimalFiber(MuscleNames_sel)/100;
    lMo_ward_SD     = getOptimalFiberSD(MuscleNames_sel)/100;
    lMo_ind         = find(isnan(lMo_ward)); % no information found
    lMo_SD_ind      = find(isnan(lMo_ward_SD)); % no information found
    
    use_ward=Misc.lMo_use_ward;
    if use_ward==0
        kT_range =[10 35];
        per_lMo=0.20;   
        per_lTs=0.20;    
        lMo_range = [(1-per_lMo)*lMo_list(ind_US) (1+per_lMo)*lMo_list(ind_US)];
        lTs_range = [(1-per_lTs)*lTs_list(ind_US) (1+per_lTs)*lTs_list(ind_US)];
    elseif use_ward==1
        kT_range =[10 35]; 
        per_SD =1;
        per_lTs=0.20;      
        lMo_range = [lMo_ward-per_SD*lMo_ward_SD lMo_ward+per_SD*lMo_ward_SD];
        lTs_range = [(1-per_lTs)*lTs_list(ind_US) (1+per_lTs)*lTs_list(ind_US)];
 
    end
    % set bounds
    opti.subject_to(kT_range(1) <kT_opt< kT_range(end));
    opti.subject_to(lMo_range(:,1) <lMo_opt< lMo_range(:,2));
    opti.subject_to(lTs_range(:,1) <lTs_opt< lTs_range(:,2));
    
    % set initial guess
    opti.set_initial(kT_opt,25);
    opti.set_initial(lMo_opt,lMo_list(ind_US));
    opti.set_initial(lTs_opt,lTs_list(ind_US));
    lMtilde_US=lMtilde(ind_US,:);
    
    % objective function
    J=J+Misc.wUS.*sumsqr(lMtilde_US-USDigitalizedInterp(2:end,:))/N/NMuscles;
    
% %     Achilles tendon - To constraint Achilles tendon to physiological range
%     [kT_abs_computed] = unNormalizedTendon(kT_opt(ind_US_AT_in),FMo_list(ind_US_AT),lTs_opt(ind_US_AT_in));
%     opti.subject_to(140 < sum(kT_abs_computed) < 171);
end

end