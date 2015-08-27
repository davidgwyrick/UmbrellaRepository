% function [sim_time, time_done, NumWorkers] = calc_sim_time (NumRuns, pCa_Range,  Pulse_Width_Range,HalfSL_Range, Rate_Range, Muscle_Type_Range, SXB2_Range, SXB3_Range, eta_range,f1_range,y0_range)
function [sim_time, time_done, NumWorkers] = calc_sim_time(NumRuns, pCa_Range, HalfSL_Range, Pulse_Width_Range, Rate_Range, Muscle_Type_Range)

time = clock;

NumWorkers = matlabpool('size');
if NumWorkers == 0
    NumWorkers = 1;
end

sim_time = 0;

if NumRuns == 0
    for ipCa = 1:length(pCa_Range)
        pCa = pCa_Range(ipCa);
        Ca = 10^(-pCa);
        maxt = SimLength(Ca);
        NSTEPS = maxt * 1000;
        TotalRuns = SimRuns(6400, 0.1, NSTEPS, maxt);
        if rem(TotalRuns,NumWorkers) > 0
            iterations = floor(TotalRuns/NumWorkers) + 1;
        else 
            iterations = floor(TotalRuns/NumWorkers);
        end 
        if pCa < 4.1
            sim_time = sim_time + (60*iterations);
        else
            sim_time = sim_time + (120*iterations);
        end
    end
else
    for ipCa = 1:length(pCa_Range)
        pCa = pCa_Range(ipCa);
        if rem(NumRuns,NumWorkers) > 0
            iterations = floor(NumRuns/NumWorkers) + 1;
        else 
            iterations = floor(NumRuns/NumWorkers);
        end
        if pCa < 4.1
            sim_time = sim_time + (60*iterations);
        else
            sim_time = sim_time + (110*iterations);
        end
    end
end

% Full_iterations = length(HalfSL_Range)*length(Rate_Range)*length(Muscle_Type_Range)*length(Pulse_Width_Range)*length(SXB2_Range)*length(SXB3_Range)*length(eta_range)*length(f1_range)*length(y0_range);
Full_iterations = length(HalfSL_Range)*length(Rate_Range)*length(Muscle_Type_Range)*length(Pulse_Width_Range);
sim_time = sim_time * Full_iterations;
time(5) = time(5) + (sim_time/60);
time(4) = time(4) + floor(time(5)/60);
time(5) = rem(time(5),60);
time_done = [num2str(time(4)), ':', num2str(time(5), '%02.0f')];

fprintf('\nTotal Estimated Time: %d minutes (around %s)\nRunning on %d cores\n',floor(sim_time/60), time_done, NumWorkers)