%% MultipleMultiple.m
% Original Author: Axel
% Updated by: David
% Last Updated: 8/19/15
%%
clc 
% Enables parallel language features like parfor
if matlabpool('size') == 0, matlabpool(4), end

TnDensitiesUniform = [];
TnDensitiesRandom = 1;
XBDensitiesUniform = [];
XBDensitiesRandom = 1;
kxscaler = 3; %cross-bridge spring scaling constant

init_params;

% SXB2_Range = [1];
% SXB3_Range = [1];

% eta_range = [0.8];
% f1_range = [0.5];
% y0_range = [35];
% SlopeL_range = [-60];%, -80, -120];
% SlopeR_range = [20];%, 24, 16];
% atan_max_range = [50];%, 80, 120];
% GausPeak_range = [2000];%, 2400, 1600];
%% Testing TF Rates (9/3/15)

Koff_range = [150];
RuOff_range = [100];
CaOff_range = [100];

%% Timer/Time Calculation
NumRuns = 12;
%pCaV = [4 4.5, 5.0, 5.5, 5.7, 5.8, 5.9, 6.0, 6.1, 6.25, 6.5, 7.0];
pCaV = [5.6];
HalfSL_Range = [1210];
Pulse_Width_Range = [160]; %number of ms of Ca2+ pulse
NumTwitch = 2; %number of twitches Max: 3 twitches
TimeBtwTwitches = [505]; %number of ms btw peak of pulses; i.e. if the pulse width range were 50ms and you
% inputed a timeBtwTwitches to be 35ms, then the start of the second twitch
% would begin as the first Ca Transient is going down
% Rate_Range = [-0.1 -0.25 -0.5 -1.0];
Rate_Range = [0];
%OutDir_Array={''};
Muscle_Type_Range = {'Soleus'};

% calc_sim_time(NumRuns, pCaV, HalfSL_Range, Rate_Range, Muscle_Type_Range, SXB2_Range, SXB3_Range, Pulse_Width_Range, eta_range,f1_range,y0_range);
calc_sim_time(NumRuns, pCaV, HalfSL_Range, Pulse_Width_Range, Rate_Range, Muscle_Type_Range);
tStart = tic;
%% Type of Ca profile for twitch analysis
%Ca_protocol = 'None';
%Ca_protocol = 'Step';
%Ca_protocol = 'Burst';
%Ca_protocol = 'Train';
Ca_protocol = 'Twitch';
%Ca_protocol = 'MultipleTwitch';

%Ca_protocol = 'LandTwitchHuman';
Ca_protocol = 'LandTwitchRat';
%Ca_protocol = 'LandTwitchMouse';
CaProfile_Scalar = 1; %0 if you do NOT want to scale Land's calcium profile transients;
%1 if you do; if you enter 1, the profile will be moved to a diastolic Ca
%concentration of pCa=8 and the peak value will be given by pCaV
%% 

% for ieta_range = 1:length(eta_range);
%     tcparam.eta = eta_range(ieta_range);
% for if1_range = 1:length(f1_range);
%     tcparam.f1 = f1_range(if1_range);
% for iy0_range = 1:length(y0_range);
%     tcparam.y0 = y0_range(iy0_range);
% for iSlopeL_range = 1:length(SlopeL_range);
%     tcparam.SlopeL = SlopeL_range(iSlopeL_range);
% for iSlopeR_range = 1:length(SlopeR_range);
%     tcparam.SlopeR = SlopeR_range(iSlopeR_range);
% for iatan_max_range = 1:length(atan_max_range);
%     tcparam.atan_max = atan_max_range(iatan_max_range);
% for iGausPeak_range = 1:length(GausPeak_range);
%     tcparam.GausPeak = GausPeak_range(iGausPeak_range);

for ikoff=1:length(Koff_range)
    Koff=Koff_range(ikoff);
for iRuoff=1:length(RuOff_range)
    RuOff=RuOff_range(iRuoff);
for iCaoff=1:length(CaOff_range)
    CaOff=CaOff_range(iCaoff);
for ipulse_width = 1:length(Pulse_Width_Range) %Loopthrough/do simulations on different Pulse widths of Ca transients
    pulse_width = Pulse_Width_Range(ipulse_width);
    %     for iSXB2 = 1:length(SXB2_Range)
    %         tcparam.SXB2=SXB2_Range(iSXB2);
    %         for iSXB3 = 1:length(SXB3_Range)
    %             tcparam.SXB3=SXB3_Range(iSXB3);
    for i = 1:length(Muscle_Type_Range) %Loopthrough/do simulations on different types of muscle
        Muscle_Type=Muscle_Type_Range{i};
        for j = 1:length(HalfSL_Range) %Loopthrough/do simulations on different sarcomere length
            StartLength = HalfSL_Range(j);
            for k = 1:length(Rate_Range) %Loopthrough/do simulations on different rates of lengthening/shortening of the sarcomere
                Rate = Rate_Range(k);
                for TnKOType=[0,1]
                    filaments.TnKOType = TnKOType;
                    %if the knockout of RUs is random
                    if TnKOType == 0, TnDensity = TnDensitiesRandom; else TnDensity = TnDensitiesUniform; end
                    for TnKO = TnDensity
                        filaments.TnFraction = TnKO;
                        for XBKOType=[0,1]
                            filaments.XBKOType = XBKOType;
                            if XBKOType == 0, XBDensity = XBDensitiesRandom; else XBDensity = XBDensitiesUniform; end
                            for XBKO = XBDensity
                                filaments.XBFraction = XBKO;
                                for kxscaler = kxscaler
                                    StiffScale.kxscaler = kxscaler;
                                    for pCa_index = 1:length(pCaV)
                                        pCa_in=pCaV(pCa_index);
                                        % pCa is a parameter by itself, i.e. no struct
                                        Timestr=Strcat(datestr(clock,'yyyy_mm_dd_HHMM'));
                                        %if strcmp(Ca_protocol,'None') == 1, OutDir = ['DataFiles' filesep 'SteadyStateCaOff=20_koff=50' filesep]; else OutDir = ['DataFiles' filesep num2str(pulse_width), 'ms', Ca_protocol filesep];end
                                        %Fix String Array later
                                        %if pCa_index ==1, OutDir_Array=OutDir; else OutDir_Array=Strcat(OutDir_Array,OutDir);end
%                                         OutDir = ['DataFiles' filesep 'TFrates_koff=',num2str(Koff), ' RuOff=',num2str(RuOff),' CaOff=',num2str(CaOff) filesep];
                                        OutDir = ['DataFiles' filesep 'CaTransient' filesep 'Scaled_',Ca_protocol,'_pCapeak=',num2str(pCa_in) filesep];
                                        mkdir(OutDir);
                                        save([OutDir filesep 'Parameters.mat']);
                                        
                                        disp( [ 'pCa = ' num2str( pCa_in )])
                                        [Steps, Stats, IndexThalf, Binder] = RunSeveral(NumRuns, DataParams, Muscle_Type, StartLength, pCa_in, StiffScale, filaments, knockout, coop, TFRateScale, tcparam, Rate, Ca_protocol, pulse_width,NumTwitch,TimeBtwTwitches,Koff,RuOff,CaOff,CaProfile_Scalar);
                                        WriteText(OutDir, pCa_in, DataParams.dt, Binder, Steps, Stats, IndexThalf);
                                        WriteTon(DataParams.dt, OutDir, pCa_in, OutDir, knockout.XB_Fraction, Stats);
%                                         clf(figure(1))
%                                         plotTSwyrick(pCaV,OutDir,Ca_protocol,pulse_width)
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    %         end
    %     end
end
end
end
end
% end
% end
% end
% end
% end
% end
% end

% Output time/force and time/FA graphs  -Axel
clf(figure(1))
plotTSwyrick(pCaV,OutDir,Ca_protocol,Pulse_Width_Range)
% Data{1,4}=[];%{pCaV,Steps,Binder,IndexThalf};
% Data=GatherData_v1(Data,OutDir_Array,pCaV)
% Process_HillCurves_v1(pCaV,OutDir,Data,3,Ca_protocol,knockout.TnFraction,knockout.XB_Fraction)

tEnd = toc(tStart);


fprintf('\nTotal Time: %d minutes and %3.2f seconds\n',floor(tEnd/60),rem(tEnd,60))