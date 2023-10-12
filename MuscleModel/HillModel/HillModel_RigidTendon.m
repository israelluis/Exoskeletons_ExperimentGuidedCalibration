function [FM, lMtilde, FMactFL, FMactFV, FMpas, cos_alpha] = HillModel_RigidTendon(a,lMT,vMT,Parameters)
% Returns muscle forces depending on muscle state assuming a rigid tendon.

% get the muscle parameters
Fmax = Parameters(1,:);
lMopt = Parameters(2,:);
lTs = Parameters(3,:);
alphaopt = Parameters(4,:);
vMmax  = Parameters(5,:);
kpe    = Parameters(6,:);
ksf    = Parameters(7,:);
ke0    = Parameters(8,:);

% Hill-type muscle model: geometric relationships
w = lMopt.*sin(alphaopt);
lM = sqrt((lMT-lTs).^2+w^2); % Rigid Tendon: lT = lTs
lMtilde = lM./lMopt;
lT = lTs*ones(size(lMtilde)); % rigid tendon
cos_alpha = (lMT-lT)./lM;
vMTtilde = vMT./lMopt;


% get the force-length- velocity characteristics
[FMpas,FMactFL,FMactFV] = getForceLengthVelocityProperties_setPassiveParam(lMtilde,vMTtilde,vMmax,kpe,ksf,ke0);

% Active muscle force
FMact = a.*FMactFL.*FMactFV;

% Muscle force
FM_norm = (FMact+FMpas);
FM = Fmax*FM_norm;

return