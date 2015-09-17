function [ROI_index] = T_OnTwitchROI(NSTEPS, Ca_protocol, pulse_width)
%Generate a matrix which contains the timesteps of where to calculate T_on
%for twitch dynamics

switch Ca_protocol
    case 'None' 
        DataSpan = 0.1; 
        EndStepStart = floor((1-DataSpan)*NSTEPS);
        ROI_index = {[EndStepStart:NSTEPS]};
    case 'Twitch'
        %Set region for analysis and returns that index
        ROI_index = {[2:pulse_width]}; % the time steps we want to look at for Twitch dynamics
    case 'LandTwitchHuman'
        ROI_index = {[2:400]}; % 400 ms is the approximate size of Land's Calcium profile
    case 'LandTwitchRat'
        ROI_index = {[2:168]};
    case 'LandTwitchMouse'
        ROI_index = {[2:268]};
end


end

