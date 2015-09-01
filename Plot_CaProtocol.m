%function Plot_Folders

%% MOD HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Calcium protocol
%protocol = 'Step';  %pCa values available: 5.5 and 8.0
%protocol = 'Burst';  %pCa values available: 4.0, 5.0, 5.5
%protocol = 'Train'; %pCa values available: 4.0, 5.0, 5.5
protocol = 'Twitch';

plot_pCa = 4.0; %See Available pCa above
Pulse_Width = 80; %30, 40, or 50, but does NOT apply to 'Step' protocol
%Trace 1
T1_SXB2 = 0.5; %0.5 or 1.0
T1_SXB3 = 0.1; %0.1 or 1.0
%Trace 2
T2_SXB2 = 1.0; %0.5 or 1.0
T2_SXB3 = 1.0; %0.1 or 1.0
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

if strcmp(protocol,'Train') == 1
    pCa_Range = [4.0, 5.0, 5.5];
elseif strcmp(protocol,'Burst') == 1
    pCa_Range = [4.0, 5.0, 5.5];
elseif strcmp(protocol,'Step') == 1
    pCa_Range = [5.5, 8.0];
elseif strcmp(protocol,'Twitch') == 1
    pCa_Range = [4.0, 5.0 5.5];
end
Pulse_Width_Range = [80]; %set to 30, 40, or 50. Or can use all 3, and then add a 4th dimension to the variables.

SXB2_Range = [0.5 1];
SXB3_Range = [0.1 1];

prcnt_Full = zeros(length(SXB2_Range),length(SXB3_Range),length(pCa_Range));  %prcnt_Full @ (SXB2 index, SXB3 index, pCa index)
ton_Full = zeros(length(SXB2_Range),length(SXB3_Range),length(pCa_Range));
ton_Back = zeros(length(SXB2_Range),length(SXB3_Range),length(pCa_Range));
ton_ms = zeros(length(SXB2_Range),length(SXB3_Range),length(pCa_Range));
Force = zeros(length(SXB2_Range),length(SXB3_Range),length(pCa_Range));
Force_std = zeros(length(SXB2_Range),length(SXB3_Range),length(pCa_Range));
TS_Time = {{},{},{}};  % TimeSeries Time @ (SXB2 index, SXB3 index, pCa index)
TS_Force = {{},{},{}}; % TimeSeries Force @ (SXB2 index, SXB3 index, pCa index)
Ca_level = {{},{},{}};


for ipulse_width = 1:length(Pulse_Width_Range)
    pulse_width = num2str(Pulse_Width_Range(ipulse_width));
    for iSXB2 = 1:length(SXB2_Range)
        SXB2 = num2str(SXB2_Range(iSXB2),'%1.2f');
        for iSXB3 = 1:length(SXB3_Range)
            SXB3 = num2str(SXB3_Range(iSXB3),'%1.2f');
            if strcmp(protocol,'Train') == 1
                foldername = ['Muscle=Soleus_Rate=0_HalfSl=1210_SXB2=',SXB2,'_SXB3=',SXB3,' ', pulse_width, 'ms Train\'];
            elseif strcmp(protocol,'Burst') == 1
                foldername = ['Muscle=Soleus_Rate=0_HalfSl=1210_SXB2=',SXB2,'_SXB3=',SXB3,' ', pulse_width, 'ms Burst\'];
            elseif strcmp(protocol,'Step') == 1
                foldername = ['Muscle=Soleus_Rate=0_HalfSl=1210_SXB2=',SXB2,'_SXB3=',SXB3,' Step\'];
            elseif strcmp(protocol,'Twitch') == 1
                foldername = ['80msTwitch\'];
            end          
            for ipCa = 1:length(pCa_Range)
                pCa = num2str(pCa_Range(ipCa),'%1.2f');
                filepath = ['J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\DataFiles\',foldername,'SS_ton_AVG_pCa_',pCa,'_ROI_A.txt'];
                %filepath = ['J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex\DataFiles\SXB2v3_CaProtocols/',foldername,'SS_ton_AVG_pCa_',pCa,'_ROI_A.txt'];
                Data = importdata(filepath);
                prcnt_Full(iSXB2,iSXB3,ipCa) = Data.data(13);
                ton_Full(iSXB2,iSXB3,ipCa) = Data.data(11);
                ton_Back(iSXB2,iSXB3,ipCa) = Data.data(15);
                ton_ms(iSXB2,iSXB3,ipCa) = Data.data(7);
                
                filepath = ['J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\DataFiles\',foldername,'SSData_pCa_',pCa,'_ROI_A.txt'];
                %filepath = ['J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex\DataFiles\SXB2v3_CaProtocols/',foldername,'SSData_pCa_',pCa,'_ROI_A.txt'];
                Data = importdata(filepath);
                Force(iSXB2,iSXB3,ipCa) = mean(Data.data(:,3));
                Force_std(iSXB2,iSXB3,ipCa) = std(Data.data(:,3));
                
                filepath = ['J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\DataFiles\',foldername,'TimeSeriesAvg_pCa_',pCa,'.txt'];
                %filepath = ['J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex\DataFiles\SXB2v3_CaProtocols/',foldername,'TimeSeriesAvg_pCa_',pCa,'.txt'];
                Data = importdata(filepath);
                TS_Time(iSXB2,iSXB3,ipCa) = {Data.data(:,1)};
                TS_Force(iSXB2,iSXB3,ipCa) = {Data.data(:,2)};
                Ca_level(iSXB2,iSXB3,ipCa) = {Data.data(:,12)};
            end
        end
    end
end

pCa_index = find(pCa_Range == plot_pCa); % index of desired pCa from pCa_Range (ie 2 = the second pCa value in the range). Optional: Create a pCa_index2 for comparison plots

iT1_SXB2 = find(SXB2_Range == T1_SXB2);
iT1_SXB3 = find(SXB3_Range == T1_SXB3);
iT2_SXB2 = find(SXB2_Range == T2_SXB2);
iT2_SXB3 = find(SXB3_Range == T2_SXB3);

Trace_1X = TS_Time{iT1_SXB2,iT1_SXB3,pCa_index};
Trace_1Y = TS_Force{iT1_SXB2,iT1_SXB3,pCa_index};
Trace_1Label = ['SXB2/3 = ', num2str(SXB2_Range(iT1_SXB2),'%1.1f'), '/', num2str(SXB3_Range(iT1_SXB3),'%1.1f')];
Trace_2X = TS_Time{iT2_SXB2,iT2_SXB3,pCa_index};
Trace_2Y = TS_Force{iT2_SXB2,iT2_SXB3,pCa_index};
Trace_2Label = ['SXB2/3 = ', num2str(SXB2_Range(iT2_SXB2),'%1.1f'), '/', num2str(SXB3_Range(iT2_SXB3),'%1.1f')];

clf(figure(1))
subplot(3,1,1)
plot(Trace_1X, Trace_1Y, ...
    Trace_2X, Trace_2Y)
legend(Trace_1Label, ...
    Trace_2Label, ...
    'location', 'Northeast')
title(['TimeSeries pCa ', num2str(pCa_Range(pCa_index),'%1.1f'), ', Pulse of ', num2str(Pulse_Width_Range), 'ms'])
xlabel('Time(s)')
ylabel('Force(pN)')

subplot(3,1,2)
[maximum max_idx] = max(Trace_1Y);
[maximum2 max_idx2] = max(Trace_2Y);
plot(Trace_1X, Trace_1Y/mean(Trace_1Y(max_idx-3:max_idx+5)), ...
    Trace_2X, Trace_2Y/mean(Trace_2Y(max_idx2-3:max_idx2+5)))
legend(Trace_1Label, ...
   Trace_2Label, ...
    'location', 'Northeast')
title(['Normalized TimeSeries pCa ', num2str(pCa_Range(pCa_index),'%1.1f'), ', Pulse of ', num2str(Pulse_Width_Range), 'ms'])
xlabel('Time(s)')
ylabel('Force(pN)')

subplot(3,1,3)
plot(Trace_1X, Ca_level{iT1_SXB2,iT1_SXB3,pCa_index})
title(['Pulse of ', num2str(Pulse_Width_Range), 'ms'])
xlabel('Time(s)')
ylabel('pCa')
set(gca, 'ydir', 'reverse')