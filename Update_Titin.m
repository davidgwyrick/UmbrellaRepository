function [M,EndVec] = Update_Titin(M,EndVec,Muscle_Type,L_TITIN,ACTNodes,NumBridges,NMYO,KMYO,LMYO,SL_EndLength,iStep)

if iStep == 1
    L_TITIN_prev = SL_EndLength(iStep)-(LMYO*(NumBridges-1));
else
L_TITIN_prev = SL_EndLength(iStep-1)-(LMYO*(NumBridges-1)); %Current SL - Thick Fil length.
end
K_TITIN_prev = TitinL2K(Muscle_Type, L_TITIN, L_TITIN_prev);


L_TITIN_Stretch = SL_EndLength(iStep)-(LMYO*(NumBridges-1)); %Current SL - Thick Fil length.
K_TITIN = TitinL2K(Muscle_Type, L_TITIN, L_TITIN_Stretch);


% K_TITIN = 0.1344;

TITIN_dif = K_TITIN - K_TITIN_prev; 

%TITIN_dif = 0;

for cntr = 1:NMYO
    EndVec(ACTNodes + (cntr-1)*NumBridges + 1) = (- KMYO * LMYO) + (K_TITIN*L_TITIN);  
    M(ACTNodes + (cntr-1)*NumBridges + 1, ACTNodes + (cntr-1)*NumBridges + 1) = M(ACTNodes + (cntr-1)*NumBridges + 1, ACTNodes + (cntr-1)*NumBridges + 1) + TITIN_dif;
end


end