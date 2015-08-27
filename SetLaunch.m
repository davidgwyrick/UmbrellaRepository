function [] = SetLaunch(Dir, SetNumberStart, kxscaler, START_LENGTH, L_TITIN)

SetNumber = SetNumberStart;

Launch(Dir, SetNumber, 0, kxscaler, 0, 1.0, START_LENGTH, L_TITIN); %TnKOType doesn't matter here

SetNumber = SetNumber + 1;

for TnKOType=[1,0] % Excel sheet has uniform first
    for TnFraction=[0.9,0.8,0.7,0.6,0.5]
        Launch(Dir, SetNumber, 0, kxscaler, TnKOType, TnFraction, START_LENGTH, L_TITIN);
        SetNumber = SetNumber + 1;
    end
end

end

