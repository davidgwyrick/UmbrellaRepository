function [fifties] = pCa50_shiv( sets )

fifties = zeros(2, length(sets));
%fifties(1, 1:length(sets)) = sets;

for iSet=1:length(sets)
    load(['DataFiles' filesep 'Set' num2str(sets(iSet)) filesep 'SimulationData.mat']);
    disp(sets(iSet));
    force50 = Process_Vert_3SimpleHillCurves_v2(pCaRange, Output_Location, OutTotal, 1, SetNumber, TnFraction, 1);
    fifties(1, iSet) = SetNumber;
    fifties(2, iSet) = force50;
end

