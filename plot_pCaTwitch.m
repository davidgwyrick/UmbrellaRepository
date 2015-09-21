%% Plotting script for Twitch Dynamics
% David Wyrick
% 9/21/15
%%
pCa_Range= [4, 4.5, 5.0, 5.25, 5.6, 5.75, 6.0];
pulse_width=[100];

%Ca_protocol = 'Twitch';
%Ca_protocol = 'MultipleTwitch';
Ca_protocol = 'LandTwitchHuman';
%Ca_protocol = 'LandTwitchRat';
%Ca_protocol = 'LandTwitchMouse';
CaProfile_Scalar = 1; %0 if you do NOT want to scale Land's calcium profile transients;
%1 if you do; if you enter 1, the profile will be moved to a diastolic Ca
%concentration of pCa=8 and the peak value will be given by pCaV
if strcmp(Ca_protocol,'LandTwitchHuman')==1,pulse_width=[400];end 
if strcmp(Ca_protocol,'LandTwitchMouse')==1 || strcmp(Ca_protocol,'LandTwitchRat')==1,pulse_width=[168];end 


TS_Time=cell(1,7);
TS_Force=cell(1,7);
TS_CaLevel=cell(1,7);
TS_XBfrac=cell(1,7); %XB Fraction Bound
TS_TF1=cell(1,7); %Actins Ca0
TS_TF2=cell(1,7); %Actins Ca1
TS_TF3=cell(1,7); %Actins Ca2
Time2Relax=cell(1,7);
T_on
index=1;
for ii=1:length(pCa_Range),
    pCa=pCa_Range(ii);
    foldername=[num2str(pulse_width),'ms_',Ca_protocol];
    filepath=['J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\DataFiles\CalciumTwitchDynamics\',foldername,'\TimeSeriesAvg_pCa_',num2str(pCa,'% 10.2f'),'.txt'];
    
    Data = importdata(filepath);
    TS_Time(index)={Data.data(:,1)};
    TS_Force(index)={Data.data(:,3)};
    TS_CaLevel(index)={Data.data(:,12)};
    TS_XBfrac(index)={Data.data(:,5)};
    TS_TF1(index)={Data.data(:,8)}; 
    TS_TF2(index)={Data.data(:,9)};
    TS_TF3(index)={Data.data(:,10)};
    
    
    filepath=['J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\DataFiles\CalciumTwitchDynamics\',foldername,'\TimeToRelax_pCa_',num2str(pCa,'% 10.2f'),'.txt'];
    Data = importdata(filepath);
    Time2Relax(index)={Data.data};
    
    filepath=['J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\DataFiles\CalciumTwitchDynamics\',foldername,'\Twitch_ton_Raw_pCa_',num2str(pCa,'% 10.2f'),'_ROI_A.txt'];
    Data = importdata(filepath);
    Time2Relax(index)={Data.data};
    index=index+1;
end
    