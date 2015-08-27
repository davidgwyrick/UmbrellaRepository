function [] = MultiLaunch(Dir, SetNumberStart, kxscaler, TnFraction, L_TITIN)

SetNumber = SetNumberStart;

for i=1:length(kxscaler)
    for j=1:length(TnFraction)
        for k=1:length(L_TITIN)
            Launch(Dir, SetNumber, 0, kxscaler(i), TnFraction(j), L_TITIN(k));
            SetNumber = SetNumber + 1;
        end
    end
end

end

