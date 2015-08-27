function [ runs ] = pCaRuns(pCa_in)

dt=(1/1)*1e-3; % simulation timestep, in seconds OH NO CONSTANT
Ca = 10^(-pCa_in); % calcium level
maxt = SimLength(Ca); % simulation length, in seconds
NSTEPS = ceil(maxt/dt); % number of simulation steps to do.
DataSpan = 0.1; % what fraction of end of sim we'll average over OH NO CONSTANT

runs = SimRuns(6400, DataSpan, NSTEPS, maxt);