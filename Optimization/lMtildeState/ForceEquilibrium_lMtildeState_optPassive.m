function [err, FT, Fpe, FMltilde, FMvtilde, cos_alpha] = ForceEquilibrium_lMtildeState_optPassive(a,...
    lMtilde,vMtilde,lM_projected,lMT,params,factor,kT,shift)
% Hill-type muscle model: equilibrium between muscle and tendon forces
% All muscle-tendon characteristics are fully described in the publication
% and its online supplement

FMo = params(:,1);
lMo = params(:,2);
lTs = params(:,3);
vMtildemax = params(:,5);
kpe = params(:,6);
so  = params(:,7);
sM  = params(:,8);

FM_factor=factor';

% Hill-type muscle model: geometric relationships
lM = lMtilde.*lMo;
lT = lMT - lM_projected;
lTtilde = lT./lTs;

% Tendon force-length characteristic
fse = (exp(kT.*(lTtilde - 0.995)))/5-0.25+shift;

% get muscle force-length characteristic
[Fpe,FMltilde,FMvtilde] = getForceLengthVelocityProperties_setPassiveParam(lMtilde,vMtilde,vMtildemax,kpe,so,sM);

% Active muscle force
Fce = a.*FMltilde.*FMvtilde;

% Muscle force
FM = Fce.*FM_factor+Fpe;

% Tendon force
FT = FMo.*fse;

% Equilibrium between muscle and tendon forces
% Fm*cos(alpha) = Ft
cos_alpha = (lMT-lT)./lM;
err =  FM.*cos_alpha-fse;

end