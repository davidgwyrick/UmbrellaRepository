sets = [604];
pCaRange=[
   4.0000
    4.5
    5.0
%    5.25
     5.5 
     5.7 
    5.8
    5.9
    6.0
    6.1
    6.25
    6.5
    7.0
    7.5
    8.0
%    9.0000
    ];
for setN = 1:length(sets)
    for pCaN = 1:length(pCaRange)
        dir = ['.' filesep 'DataFiles' filesep 'Set' num2str(sets(setN)) filesep];
        try
            ton_pCa_Data(dir,pCaRange(pCaN),dir,sets(setN));
        catch
            disp(['Ignoring set #' num2str(sets(setN)) ' at pCa ' num2str(pCaRange(pCaN))]);
        end
    end
end