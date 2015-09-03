function [Steps, Stats, IndexThalf, Binder] = RunSeveral(RequestedRuns, DataParams, Muscle_Type, StartLength, pCa, StiffScale, filaments, knockout, coop, TFRateScale, tcparam, Rate, Ca_protocol, pulse_width, NumTwitch, TimeBtwTwitches,Koff,RuOff,CaOff)

%% Calcium Level in Molar
Ca = 10^(-pCa); 

%% Data/sim level stuff
DataSpan = DataParams.DataSpan;
dt=DataParams.dt; % simulation timestep, in seconds
maxt = SimLength(Ca); % simulation length, in seconds
NSTEPS = ceil(maxt/dt); % number of simulation steps to do.

%% SETUP Length change index and Regions of interest (ROIs) for analysis 
%Rate = 0; %Length change in (ML/s). (-) for shortening;(+) for lengthing.
%SS_Steps = NSTEPS;  %Hold on to the number of steps before any length change occurs to be used in SS analysis (just returned by the function)
[SL_EndLength, NSTEPS, ROI_index] = LengthChange(StartLength, NSTEPS, Rate, DataParams); %Returns index of handle position for each step, as well as # of steps to add; Assumes thin filiment length of 1119nM.

%% SETUP Ca level change index 
[CaChange_index] = CaChange(Ca_protocol, pCa, NSTEPS, pulse_width, NumTwitch, TimeBtwTwitches);

%% Calculates the total number of runs 
if (RequestedRuns == 0),TotalRuns = SimRuns(6400, DataSpan, NSTEPS, maxt); else TotalRuns = RequestedRuns; end

%% Knockout!

TnKOType = knockout.TnKOType; % 0 for random, 1 for uniform
XBKOType = knockout.XBKOType; % ditto

TnFraction=knockout.TnFraction; % Functional Tn density
XB_Fraction=knockout.XB_Fraction; % Functional XB density

%% Filament and stiffness info

kxscaler = StiffScale.kxscaler; %Cross-Bridge Spring scalar
Tm_Type = filaments.Tm_Type; %
LACT = filaments.LACT;
LMYO = filaments.LMYO;
NACT = filaments.NACT; % N IS 'NUMBER OF' filaments actin
NMYO = filaments.NMYO; % N IS 'NUMBER OF' filaments myosin
NumBridges = filaments.NumBridges; % bridges per filament + handel
NumSites = filaments.NumSites; % sites per filament + handel
N_Thick_Start = filaments.N_Thick_Start; %% 3-start thick filaments (3 heads per crown)
L_TITIN = filaments.L_TITIN;
L_TITIN_Stretch = StartLength - (LMYO*(NumBridges-1));
NTn = filaments.NTn;

ACTNodes=NACT*NumSites;
MYONodes=NMYO*NumBridges;
TOTNodes=ACTNodes+MYONodes;

KACT = StiffScale.kascaler*3*filaments.KACT_0;
KMYO = StiffScale.kmscaler*3*filaments.KMYO_0;

%K_TITIN = 0; 
K_TITIN = TitinL2K(Muscle_Type, L_TITIN, L_TITIN_Stretch);

%% Cooperativity
CoopPass=zeros(13, 1);
CoopPass(1, 1)=coop.TF3TF12;
CoopPass(2, 1)=coop.XB2TF12;
CoopPass(3, 1)=coop.XB3TF12;
CoopPass(4, 1)=coop.XB2TF23;
CoopPass(5, 1)=coop.XB3TF23;
CoopPass(6, 1)=coop.TF3TF23;
CoopPass(7, 1)=coop.XB2TF21;
CoopPass(8, 1)=coop.XB3TF21;
CoopPass(9, 1)=coop.TF3TF21;
CoopPass(10, 1)=coop.XB2TF32;
CoopPass(11, 1)=coop.XB3TF32;
CoopPass(12, 1)=coop.TF3TF32;
CoopPass(13, 1)=coop.Coop_Type;

%% thermochemical stuff
thermochem.SXB1=tcparam.SXB1;
thermochem.SXB2=tcparam.SXB2;
thermochem.SXB3=tcparam.SXB3;

thermochem.eta = tcparam.eta; % The mechanical efficiency fraction for the total amount of work performed.
thermochem.f1 = tcparam.f1; % The fraction drop from the first energy to the second minimum
%thermo.f1 = 0.28*(1/XB_Fraction);
% f2 = eta - f1; % To complete the balance of drop from second E min. to the 3rd E min.

y0=tcparam.y0; % Offset
xL=tcparam.xL; % xposition left of 0;
xR=tcparam.xR; % xposition right of 0;
yL=tcparam.SlopeL*xL + y0; % yposition left of 0;
yR=tcparam.SlopeR*xR + y0; % yposition right of 0;
m=(-y0 + yR - (0.5*sqrt( ((xL*(y0-yR)-xR*(y0-yL))/xL)^2 )))/xR;
A=((xL*(y0-yR)-xR*(y0-yL))^2)/(4*(xL^2)*(xR^2));

thermochem.AA=1;   %% k12 offset, not sure if should be parametrized
thermochem.BB=A;   %% k20 ???
thermochem.CC=m;   %% k20 ???
thermochem.DD=y0;  %% k20 offset--right?

thermochem.atan_max = tcparam.atan_max;
thermochem.GausPeak = tcparam.GausPeak;

J2pNnm=1e21;
kCal2Joule=4.1868e3;
N_Avo=6.022e23;
% Gnot = 13*4.1868e3*1e21/6.022e23; % 13 kCal/mol (E of hyrolysis of ATP) * kCal2Joule * J2pNnm / N_Avagadro's (molecules/mol)
% % Gnot then becomes free energy of hydrolysis in pN*nm from the above line
Gnot=7.8*kCal2Joule; % Gnot of ATP hydrolysis in Joule/mol
% Gnot=13*kCal2Joule; % Gnot of ATP hydrolysis in Joule/mol
T=273.15; % Temperature in K
J2RT=1/(8.314*T); %To convert J/mol to units of RT [R=8.314 J/(mol K) ]
GnotRT=Gnot*J2RT; %Gnot in units of RT ~ 21
thermochem.dGtot=-GnotRT-log(tcparam.ATP/(tcparam.ADP*tcparam.phos)); %Total free energy drop from the
% hydrolysis and concentration changes in units of RT
BridgeNRG=abs(thermochem.eta*thermochem.dGtot); % The total amount of energy allowed for the XB displacement
% given the efficiency term on the total free energy available.  In units
% of RT.

RTnm2 = 8.314*T*J2pNnm/N_Avo; %  to convert RT/nm^2 into pN/nm, such that 1 RT/nm^2 =~ 4 pN/nm @ 310K
thermochem.kRT = kxscaler/RTnm2; % cross brige spring constant in RT/nm^2
thermochem.reach = sqrt(BridgeNRG./thermochem.kRT); % The cross-brige reach, in nm

%% Construct some big stuff
Tn = MakeTn(NTn, NACT, NumSites);

Mcore = MakeMcore(NACT, KACT, ACTNodes, NMYO, KMYO, NumSites, NumBridges, K_TITIN, TOTNodes);
EndVeccore = MakeEndVeccore(NACT, NumSites, NMYO, NumBridges, KACT, LACT, KMYO, LMYO, StartLength, K_TITIN, L_TITIN, ACTNodes, TOTNodes);

HARDY = MakeHARDY(); % more moved out for space than anything.
LAUREL = MakeLAUREL(filaments.Angle_Crowns, NumBridges, filaments.Angle_Thick_Start, MYONodes, N_Thick_Start, NMYO);

[TFRatePass(1,1), TFRatePass(2,1), TFRatePass(3,1), TFRatePass(4,1), TFRatePass(5,1), TFRatePass(6,1)] = ScaleThinFilRates(TFRateScale.ScaleK1, TFRateScale.ScaleK2, TFRateScale.ScaleK3,Koff,RuOff,CaOff);

%% XB binding stuff
Big_Binder = {};

%% Data per time step, averaged across runs
MFvec = zeros(1,NSTEPS); %Myosin Fraction vector?
AFvec = zeros(1,NSTEPS); %Actin Fraction vector?   
FractXB1 = zeros(1,NSTEPS); %Fraction of myosin binding sites in unbound state
FractXB2 = zeros(1,NSTEPS); %fraction of myosin binding sites in bound pre-power state
FractCa0 = zeros(1,NSTEPS); %Fraction of ...
FractCa1 = zeros(1,NSTEPS);
FractCa2 = zeros(1,NSTEPS);
ATPuse = zeros(1,NSTEPS);



%% Means and Vars
Stats.MFMean = zeros(1,TotalRuns,(length(ROI_index)));
Stats.AFMean = zeros(1,TotalRuns,(length(ROI_index)));
Stats.FractXB1Mean = zeros(1,TotalRuns,(length(ROI_index)));
Stats.FractXB2Mean = zeros(1,TotalRuns,(length(ROI_index)));
Stats.FractCa0Mean = zeros(1,TotalRuns,(length(ROI_index)));
Stats.FractCa1Mean = zeros(1,TotalRuns,(length(ROI_index)));
Stats.FractCa2Mean = zeros(1,TotalRuns,(length(ROI_index)));
Stats.ATPuseMean = zeros(1,TotalRuns,(length(ROI_index)));
Stats.FractBoundMean = zeros(1,TotalRuns,(length(ROI_index)));
Stats.MFVar = zeros(1,TotalRuns,(length(ROI_index)));
Stats.AFVar = zeros(1,TotalRuns,(length(ROI_index)));
Stats.FractXB1Var = zeros(1,TotalRuns,(length(ROI_index)));
Stats.FractXB2Var = zeros(1,TotalRuns,(length(ROI_index)));
Stats.FractCa0Var = zeros(1,TotalRuns,(length(ROI_index)));
Stats.FractCa1Var = zeros(1,TotalRuns,(length(ROI_index)));
Stats.FractCa2Var = zeros(1,TotalRuns,(length(ROI_index)));
Stats.ATPuseVar = zeros(1,TotalRuns,(length(ROI_index)));
Stats.FractBoundVar = zeros(1,TotalRuns,(length(ROI_index)));

%% Time to reach half the mean (I think)
IndexThalfMFvec = zeros(1,TotalRuns);
IndexThalfAFvec = zeros(1,TotalRuns);
IndexThalfFractXB1 = zeros(1,TotalRuns);
IndexThalfFractXB2 = zeros(1,TotalRuns);
IndexThalfFractCa0 = zeros(1,TotalRuns);
IndexThalfFractCa1 = zeros(1,TotalRuns);
IndexThalfFractCa2 = zeros(1,TotalRuns);
IndexThalfATPuse = zeros(1,TotalRuns);
IndexThalfFractBound = zeros(1,TotalRuns);


FinalFractCa0 = {}; 
FinalFractCa1 = {}; 
FinalFractCa2 = {}; 
FinalFractXB1 = {}; 
FinalFractXB2 = {}; 
FinalMFvec = {}; 
FinalAFvec = {}; 
FinalATPuse = {}; 
%%

disp(['Starting... (' num2str( TotalRuns ) ' runs)'])

% for iTotRun=1:TotalRuns
parfor iTotRun=1:TotalRuns
    tRunStart = tic;
       
    [Binder, TempFractCa0, TempFractCa1, TempFractCa2, TempFractXB1, TempFractXB2, TempMFvec, TempAFvec, TempATPuse] = OneRun(Muscle_Type, Tm_Type, dt, NSTEPS, L_TITIN, LACT, KACT, NACT, LMYO, KMYO, NMYO, NTn, Tn, TnFraction, HARDY, LAUREL, XB_Fraction, ACTNodes, MYONodes, NumSites, NumBridges, Mcore, EndVeccore, N_Thick_Start, CoopPass, TFRatePass, Ca, 1, kxscaler, thermochem, TnKOType, XBKOType, SL_EndLength, CaChange_index, ROI_index);

    Big_Binder = [Big_Binder, {Binder}];

    MFvec=MFvec+TempMFvec/TotalRuns;
    AFvec=AFvec+TempAFvec/TotalRuns;
    FractXB1=FractXB1+TempFractXB1/TotalRuns;
    FractXB2=FractXB2+TempFractXB2/TotalRuns;
    FractCa0=FractCa0+TempFractCa0/TotalRuns;
    FractCa1=FractCa1+TempFractCa1/TotalRuns;
    FractCa2=FractCa2+TempFractCa2/TotalRuns;
    ATPuse=ATPuse+TempATPuse/TotalRuns;
    
    FinalFractCa0 = [FinalFractCa0, {TempFractCa0}];
    FinalFractCa1 = [FinalFractCa1, {TempFractCa1}];
    FinalFractCa2 = [FinalFractCa2, {TempFractCa2}];
    FinalFractXB1 = [FinalFractXB1, {TempFractXB1}];
    FinalFractXB2 = [FinalFractXB2, {TempFractXB2}];
    FinalATPuse = [FinalATPuse, {TempATPuse}];
    FinalMFvec = [FinalMFvec, {TempMFvec}];
    FinalAFvec = [FinalAFvec, {TempAFvec}];

    tRunEnd = toc(tRunStart);
    fprintf('Finished run %d of %d in %d minutes and %3.2f seconds\n',iTotRun,TotalRuns,floor(tRunEnd/60),rem(tRunEnd,60))
end

Binder = Big_Binder;

OutFractCa0 = FinalFractCa0; 
OutFractCa1 = FinalFractCa1 ; 
OutFractCa2 = FinalFractCa2 ; 
OutFractXB1 = FinalFractXB1; 
OutFractXB2 = FinalFractXB2; 
OutATPuse = FinalATPuse; 
OutMFvec = FinalMFvec; 
OutAFvec = FinalAFvec; 

for j = 1:TotalRuns
    for i = 1:length(ROI_index)
        Stats.MFMean(1,j,i) = mean(OutMFvec{j}(ROI_index{i}));
        Stats.AFMean(1,j,i) = mean(OutAFvec{j}(ROI_index{i}));
        Stats.FractXB1Mean(1,j,i) = mean(OutFractXB1{j}(ROI_index{i}));
        Stats.FractXB2Mean(1,j,i) = mean(OutFractXB2{j}(ROI_index{i}));
        Stats.FractCa0Mean(1,j,i) = mean(OutFractCa0{j}(ROI_index{i}));
        Stats.FractCa1Mean(1,j,i) = mean(OutFractCa1{j}(ROI_index{i}));
        Stats.FractCa2Mean(1,j,i) = mean(OutFractCa2{j}(ROI_index{i}));
        Stats.ATPuseMean(1,j,i) = mean(OutATPuse{j}(ROI_index{i}));
        Stats.FractBoundMean(1,j,i) = Stats.FractXB1Mean(1,j,i) + Stats.FractXB2Mean(1,j,i); % sums: convenient
        
        Stats.MFVar(1,j,i)=var(OutMFvec{j}(ROI_index{i}:ROI_index{i}(length(ROI_index{i}))));
        Stats.AFVar(1,j,i)=var(OutAFvec{j}(ROI_index{i}:ROI_index{i}(length(ROI_index{i}))));
        Stats.FractXB1Var(1,j,i)=var(OutFractXB1{j}(ROI_index{i}:ROI_index{i}(length(ROI_index{i}))));
        Stats.FractXB2Var(1,j,i)=var(OutFractXB2{j}(ROI_index{i}:ROI_index{i}(length(ROI_index{i}))));
        Stats.FractCa0Var(1,j,i)=var(OutFractCa0{j}(ROI_index{i}:ROI_index{i}(length(ROI_index{i}))));
        Stats.FractCa1Var(1,j,i)=var(OutFractCa1{j}(ROI_index{i}:ROI_index{i}(length(ROI_index{i}))));
        Stats.FractCa2Var(1,j,i)=var(OutFractCa2{j}(ROI_index{i}:ROI_index{i}(length(ROI_index{i}))));
        Stats.ATPuseVar(1,j,i)=var(OutATPuse{j}(ROI_index{i}:ROI_index{i}(length(ROI_index{i}))));
        Stats.FractBoundVar(1,j,i)=var(OutFractXB1{j}(ROI_index{i}) + OutFractXB2{j}(ROI_index{i})); % these are correlated probably, so we can't just sum variances
    end
end


% Also calculate the time to reach 1/2 of the SS value for each of the
% variables.
% Load them with time value, closest to half of the max avg value
% This is irritatingly verbose but I don't know how to fix it.   %%DO WE WANT TO REPEAT FOR MULTIPLE ROIs?????? -Axel
for i = 1:TotalRuns
    [~,Index] = min(abs(OutMFvec{i}-(0.5*Stats.MFMean(1,i,1))));
    IndexThalfMFvec(i)=dt*Index;
    [~,Index] = min(abs(OutAFvec{i}-(0.5*Stats.AFMean(1,i,1))));
    IndexThalfAFvec(i)=dt*Index;
    [~,Index] = min(abs(OutFractXB1{i}-(0.5*Stats.FractXB1Mean(1,i,1))));
    IndexThalfFractXB1(i)=dt*Index;
    [~,Index] = min(abs(OutFractXB2{i}-(0.5*Stats.FractXB2Mean(1,i,1))));
    IndexThalfFractXB2(i)=dt*Index;
    [~,Index] = min(abs(OutFractCa0{i}-(0.5*Stats.FractCa0Mean(1,i,1))));
    IndexThalfFractCa0(i)=dt*Index;
    [~,Index] = min(abs(OutFractCa1{i}-(0.5*Stats.FractCa1Mean(1,i,1))));
    IndexThalfFractCa1(i)=dt*Index;
    [~,Index] = min(abs(OutFractCa2{i}-(0.5*Stats.FractCa2Mean(1,i,1))));
    IndexThalfFractCa2(i)=dt*Index;
    [~,Index] = min(abs(OutATPuse{i}-(0.5*Stats.ATPuseMean(1,i,1))));
    IndexThalfATPuse(i)=dt*Index;
    [~,Index] = min(abs((OutFractXB1{i}) + OutFractXB2{i})- (0.5*Stats.FractBoundMean(1,i,1)));
    IndexThalfFractBound(i)=dt*Index;
end


Steps = [MFvec; AFvec; FractXB1; FractXB2; FractCa0; FractCa1; FractCa2; ATPuse; SL_EndLength; CaChange_index'];
IndexThalf = [IndexThalfMFvec; IndexThalfAFvec; IndexThalfFractXB1; IndexThalfFractXB2; IndexThalfFractCa0; IndexThalfFractCa1; IndexThalfFractCa2; IndexThalfATPuse; IndexThalfFractBound];

end


    
    
    
    
    
    
    

