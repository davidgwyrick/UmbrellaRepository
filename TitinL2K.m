function [K_TITIN] = TitinL2K(Muscle_Type,L_TITIN, L_TITIN_Stretch)
%Passive titin force can be represented as:
%      Fpas(Xhs) = sigma*exp((Xhs-Xoffset)/L)   %From Campbell 2009
%
%   Force for cross sectional area


%%Static Model
%K_TITIN = 0.1344 pN nm^-1
%L_TITIN = 247nm for 1/2SL of 1150nM
%
% Resulting in passive force of 25 pN (4.5 kPa for 5600nm^2 cross sectional
% area)

%%Bert Values
%Remember 24 titin molecules per simulated half-sarcomere
%   (6 titin per thick-filament (Liversage,2001) x 4 thick filaments)
%
%Say passive tension = 0 when SL = 1.9 (1/2SL = 950nm)
%Since LMYO = 860 and L_TITIN = (1/2SL - LMYO)
%Passive force = 0 when L_TITIN = 950-860 = 90nm

%% Values Taken From Titin_Curve_Fitting.m

Xhs = L_TITIN_Stretch;
Rest_length = L_TITIN; %Rest Length of titin in nm

switch Muscle_Type 
    
    case 'Soleus' %Soleus
        sigma = 0.4037;
        Xoffset = Rest_length;
        L = 200.2494;
        
    case 'Psoas' %Psoas
        sigma =  0.0776;
        Xoffset = Rest_length;
        L = 90.5242;       
        
    case 'N2B' %N2B
        sigma = 0.3713;
        Xoffset = Rest_length;
        L = 63.8005;
end
Fpas = (sigma*exp((Xhs-Xoffset)/L));
K_TITIN = Fpas/(L_TITIN_Stretch-Rest_length);

end