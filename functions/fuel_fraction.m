function [Mff_stat,Mff_breg,MLNDoMTO,WTOoWcr] = fuel_fraction(TLARS)
%UNTITLED2 Summary of this funcjion goes here
%   1 - engine start and warm-up
%   2 - taxi
%   3 - take-off
%   4 - climb (statistical)
%   5 - cruise (for given range)
%   6 - loiter (in terms of endurance)
%   7 - descend
%   8 - Landing, Taxi and Shutdown
%   9 - Alternate (for a given range)
% Definizione del Profilo di Volo

Wfrac(1) = 0.99;
Wfrac(2) = 0.99;
Wfrac(3)= 0.995;
Wfrac(4) = 0.980;
Wfrac(7) = 0.99;
Wfrac(8) = 0.992;

% Breguet
Wfrac(5) = exp( - TLARS.cruise.R*TLARS.cruise.cj/( TLARS.cruise.V * TLARS.cruise.E ) ); % Cruise
Wfrac(6) = exp( - TLARS.loiter.cj*TLARS.loiter.End/TLARS.loiter.E ); % Loiter
Wfrac(9) = exp( - TLARS.alter.R*TLARS.alter.cj/( TLARS.alter.V * TLARS.alter.E ) ); % Alternate
n_phases = 9;

Mff_stat = 1;
for i = 1:n_phases
    Mff_stat = Mff_stat*Wfrac(i);
end

TLARS.climb.End = TLARS.cruise.h/TLARS.climb.RoC; TLARS.climb.R = TLARS.climb.End*sqrt( TLARS.climb.V^2 - TLARS.climb.RoC^2 );
Wfrac(4) = exp( - TLARS.climb.cj*TLARS.climb.End/TLARS.climb.E ); % Climb

Mff_breg = 1;
for i = 1:n_phases
    Mff_breg = Mff_breg*Wfrac(i);
end
MLNDoMTO = Mff_breg/Wfrac(8);
WTOoWcr(1) = Wfrac(1)*Wfrac(2)*Wfrac(3)*Wfrac(4);
WTOoWcr(2) = Wfrac(1)*Wfrac(2)*Wfrac(3)*Wfrac(4)*Wfrac(5);
end

