function [TLARS] = read_TLARs(tlars_file_path)
%read_TLARs function that reads TLARs and assumptions for a design airplane
%from a TXT file
%   tlars_file_path: string containing TLARs txt file path
%   TLARS: struct containing the TLARs
%
%   02/05/25 v. 0.1e

f_id = fopen(tlars_file_path,'r');
% Grafica di Lettura
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Reading TLARs ...' )
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
try
    temp = fgetl(f_id); tag = 'TLARS';
    if ~strcmp(temp,tag)
        error(' The file opened does not contain TLARs')
    end
    % General TLARs
    temp = fgetl(f_id); tag = 'General';
    if ~strcmp(temp,tag)
        error(' Expected General TLARS')
    end
    
    TLARS.npax  = fscanf(f_id,'%f '); temp = fgetl(f_id);
    TLARS.ncrew   = fscanf(f_id,'%f ');
    if isempty( TLARS.ncrew ) || TLARS.ncrew < ceil( TLARS.npax/50 )
        TLARS.ncrew = ceil( TLARS.npax/50 );
    end
    temp = fgetl(f_id);
    TLARS.npil    = fscanf(f_id,'%f '); temp = fgetl(f_id);
    
    temp = fgetl(f_id); tag = 'Aerodynamics';
    if ~strcmp(temp,tag)
        error(' Expected Aerodynamics TLARS')
    end
    TLARS.AR  = fscanf(f_id,'%f '); temp = fgetl(f_id);
    TLARS.e  = fscanf(f_id,'%f '); temp = fgetl(f_id); 
    TLARS.dCD0_f_TO  = fscanf(f_id,'%f '); temp = fgetl(f_id);
    TLARS.dCD0_f_LND  = fscanf(f_id,'%f '); temp = fgetl(f_id);
    temp  = fscanf(f_id,'%f '); 
    if isempty( temp )
        temp = ( TLARS.dCD0_f_TO+TLARS.dCD0_f_LND )*0.5;
    end
    TLARS.dCD0_f_App = temp; temp = fgetl(f_id);
    TLARS.dCD0_lgs  = fscanf(f_id,'%f '); temp = fgetl(f_id);
    TLARS.de_TO  = -fscanf(f_id,'%f '); temp = fgetl(f_id);
    TLARS.de_LND = -fscanf(f_id,'%f '); temp = fgetl(f_id);
    
    temp = fgetl(f_id); tag = 'Propulsion';
    if ~strcmp(temp,tag)
        error(' Expected Propulsive TLARS')
    end
    TLARS.T0oTmc  = fscanf(f_id,'%f '); temp = fgetl(f_id);
    TLARS.nengine = fscanf(f_id,'%f '); temp = fgetl(f_id);
    
    % Take-Off
    temp = fgetl(f_id); tag = 'Take-Off';
    if ~strcmp(temp,tag)
        error(' Expected Take-Off Requirements')
    end
    temp = fscanf(f_id,'%f '); TLARS.TO.fieldmax = convlength(temp,'ft','m');
    temp = fgetl(f_id);
    
    % Climb Specs
    temp = fgetl(f_id); tag = 'Climb';
    if ~strcmp(temp,tag)
        error(' Expected Climb Requirements')
    end
    temp = fscanf(f_id,'%f '); TLARS.climb.RoC = convvel(temp,'ft/min','m/s');
    temp = fgetl(f_id);
    temp = fscanf(f_id,'%f '); TLARS.climb.V = convvel(temp,'kts','m/s');
    temp = fgetl(f_id);
    TLARS.climb.E   = fscanf(f_id,'%f '); temp = fgetl(f_id);
    temp  = fscanf(f_id,'%f '); TLARS.climb.cj = temp/3600;
    temp= fgetl(f_id);
    
    % Cruise Specs.
    temp = fgetl(f_id); tag = 'Cruise';
    if ~strcmp(temp,tag)
        error(' Expected Cruise Requirements')
    end
    temp = fscanf(f_id,'%f '); TLARS.cruise.h = convlength(temp,'ft','m');
    temp = fgetl(f_id);
    temp = fscanf(f_id,'%f '); TLARS.cruise.R = convlength(temp,'naut mi','m');
    temp = fgetl(f_id);
    [T, a_sound, P, rho] = atmosisa(TLARS.cruise.h);
    TLARS.cruise.M = fscanf(f_id,'%f '); temp = fgetl(f_id);
    TLARS.cruise.V = TLARS.cruise.M*a_sound;
    TLARS.cruise.E = fscanf(f_id,'%f '); temp = fgetl(f_id);
    temp  = fscanf(f_id,'%f '); TLARS.cruise.cj = temp/3600;
    temp = fgetl(f_id);
    
    % Loiter Specs.
    temp = fgetl(f_id); tag = 'Loiter';
    if ~strcmp(temp,tag)
        warning('Loiter not Defined')
        TLARS.loiter.End = 0;
    else
        TLARS.loiter.E = fscanf(f_id,'%f '); temp = fgetl(f_id);
        temp  = fscanf(f_id,'%f '); TLARS.loiter.cj = temp/3600;
        temp = fgetl(f_id);
        temp  = fscanf(f_id,'%f '); TLARS.loiter.End = temp*60;
        temp= fgetl(f_id);
    end
    % Alternate Specs
    temp = fgetl(f_id); tag = 'Alternate';
    if ~strcmp(temp,tag)
        warning('Alternate not Defined')
        TLARS.alter.R = 0;
        
        TLARS.alter.h = 1; TLARS.alter.M = 1; % continua
        
        
    else
        temp = fscanf(f_id,'%f '); TLARS.alter.h = convlength(temp,'ft','m');
        temp = fgetl(f_id);
        temp = fscanf(f_id,'%f '); TLARS.alter.R = convlength(temp,'naut mi','m');
        temp = fgetl(f_id);
        [T, a_sound, P, rho] = atmosisa(TLARS.alter.h);
        TLARS.alter.M = fscanf(f_id,'%f '); temp = fgetl(f_id);
        TLARS.alter.V = TLARS.alter.M*a_sound;
        TLARS.alter.E = fscanf(f_id,'%f '); temp = fgetl(f_id);
        temp  = fscanf(f_id,'%f '); TLARS.alter.cj = temp/3600;
        temp = fgetl(f_id);
        
    end
    
    % Landing Specs
    temp = fgetl(f_id); tag = 'Landing';
    if ~strcmp(temp,tag)
        error('Expected Landing TLARS')
    end
    temp = fscanf(f_id,'%f '); TLARS.LND.SGmax = convlength(temp,'ft','m');
    temp = fgetl(f_id);
    
    fclose(f_id);
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
catch ERR
    fclose('all');
    error( ['An error occurred while reading TLARs, cannot proceed ...',ERR.message] );
end


end

