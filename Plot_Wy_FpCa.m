
pCa_Range= [4, 4.5, 5.0, 5.25, 5.6, 5.75, 6.0];
pulse_width=[100];
Ca_protocol='LandHumanTwitch';
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


% prcnt_Full = zeros(1,length(pCa_Range));
% stdv_prcnt_Full = zeros(1,length(pCa_Range));
% ton_Full = zeros(1,length(pCa_Range));
% stdv_ton_Full = zeros(1,length(pCa_Range));
% ton_Back = zeros(1,length(pCa_Range));
% stdv_ton_Back = zeros(1,length(pCa_Range));
% ton_ms = zeros(1,length(pCa_Range));
% stdv_ton_ms=zeros(1,length(pCa_Range));
% Force = zeros(1,length(pCa_Range));
% Force_std = zeros(1,length(pCa_Range));
prcnt_Full = {{},{},{}};
prcnt_Back = {{},{},{}};
ton_Full = {{},{},{}};
ton_Back = {{},{},{}};
ton_ms = {{},{},{}};
XBF_Raw = {{},{},{}};
XBF_Full = {{},{},{}};
XBF_Back = {{},{},{}};
XBF_EnRaw = {{},{},{}};
XBF_EnFull = {{},{},{}};
XBF_EnBack = {{},{},{}};
Force = {{},{},{}};

TS_Time = {{},{},{}};  % TimeSeries Time @ (SXB2 index, SXB3 index, pCa index)
TS_Force = {{},{},{}}; % TimeSeries Force @ (SXB2 index, SXB3 index, pCa index)
Ca_level = {{},{},{}};

%addpath('J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex')
% foldername = ['Vary_pCa_10msTwitch\'];
foldername = ['80msTwitch\'];
%foldername = ['400msTwitch\'];


%foldername = ['SSdataTest\'];
for ipCa = 1:length(pCa_Range)
    pCa = num2str(pCa_Range(ipCa),'%1.2f');
    %filepath = ['J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex\DataFiles/',foldername,'SS_ton_AVG_pCa_',pCa,'_ROI_A.txt'];
    filepath = ['J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\DataFiles\',foldername,'SS_ton_Raw_pCa_',pCa,'_ROI_A.txt'];
    Data = importdata(filepath);
    
    prcnt_Full(ipCa) = {Data.data(:,7)};
    prcnt_Back(ipCa) = {Data.data(:,9)};
    ton_Full(ipCa) = {Data.data(:,6)};
    ton_Back(ipCa) = {Data.data(:,8)};
    ton_ms(ipCa) = {Data.data(:,4)};
    XBF_Raw(ipCa) = {Data.data(:,10)};
    XBF_Full(ipCa) = {Data.data(:,11)};
    XBF_Back(ipCa) = {Data.data(:,12)};
    XBF_EnRaw(ipCa) = {Data.data(:,13)};
    XBF_EnFull(ipCa) = {Data.data(:,14)};
    XBF_EnBack(ipCa) = {Data.data(:,15)};

    filepath = ['J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\DataFiles\',foldername,'SSData_pCa_',pCa,'_ROI_A.txt'];
    Data = importdata(filepath);
%     Force(1,ipCa) = mean(Data.data(:,3));
%     Force_std(1,ipCa) = std(Data.data(:,3));
    
    filepath = ['J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\DataFiles\',foldername,'TimeSeriesAvg_pCa_',pCa,'.txt'];
    Data = importdata(filepath);
    TS_Time(ipCa) = {Data.data(:,1)};
    TS_Force(ipCa) = {Data.data(:,2)};
    Ca_level(ipCa) = {Data.data(:,12)};
end

%% Timeseries Curves
Trace_1X=TS_Time{1};
Trace_1Y=TS_Force{1};
Trace_1Label = ['pCa=',num2str(pCa_Range(1))];
Trace_2X=TS_Time{2};
Trace_2Y=TS_Force{2};
Trace_2Label = ['pCa=',num2str(pCa_Range(2))];
Trace_3X=TS_Time{3};
Trace_3Y=TS_Force{3};
Trace_3Label = ['pCa=',num2str(pCa_Range(3))];

%Force/time curves overplot
clf(figure(2))
subplot(3,1,1);
hold on;
plot(Trace_1X,Trace_1Y,'k-'),hold on
plot(Trace_2X,Trace_2Y,'r-'),hold on
plot(Trace_3X,Trace_3Y,'g-'),hold on
xlabel('Time(s)')
ylabel('Force (pN)')
title('10ms Twitch at different pCa')
legend(Trace_1Label,Trace_2Label,Trace_3Label)
%axis([-.1 1.5 0 800],'on')
axis([-.1 2.0 0 1000],'on')
set(gca,'YTick',0:200:800)
%% Finding Time to 50% and 90% Relaxation

RT90_found=0;
RT50_found=0;
RT90_ind=0;
RT50_ind=0;
for ii=1:length(Trace_1Y),
    [mm mm_index] =max(Trace_1Y);
    if ii > mm_index && (Trace_1Y(ii)/mm)< .50 && RT50_found == 0; 
        RT50_ind=ii;
        RT50_found=1;
    end
    if ii > mm_index && (Trace_1Y(ii)/mm)< .10 && RT90_found == 0; 
        RT90_ind=ii;
        RT90_found=1;
    end
    
end
disp(['Time to 50% Relaxation: ' num2str(Trace_1X(RT50_ind)) 's for a pCa of ' num2str(pCa_Range(1)) ' at ' num2str(Pulse_Width) 'ms Twitch'])
disp(['Time to 90% Relaxation: ' num2str(Trace_1X(RT90_ind)) 's for a pCa of ' num2str(pCa_Range(1)) ' at ' num2str(Pulse_Width) 'ms Twitch'])
disp(' ')
RT90_found=0;
RT50_found=0;
RT90_ind2=0;
RT50_ind2=0;
for ii=1:length(Trace_2Y),
    [mm2 mm_index2] =max(Trace_2Y);
    if ii > mm_index2 && (Trace_2Y(ii)/mm2)< .50 && RT50_found == 0; 
        RT50_ind2=ii;
        RT50_found=1;
    end
    if ii > mm_index2 && (Trace_2Y(ii)/mm2)< .10 && RT90_found == 0; 
        RT90_ind2=ii;
        RT90_found=1;
    end
    
end
disp(['Time to 50% Relaxation: ' num2str(Trace_2X(RT50_ind2)) 's for a pCa of ' num2str(pCa_Range(2)) ' at ' num2str(Pulse_Width) 'ms Twitch'])
disp(['Time to 90% Relaxation: ' num2str(Trace_2X(RT90_ind2)) 's for a pCa of ' num2str(pCa_Range(2)) ' at ' num2str(Pulse_Width) 'ms Twitch'])
disp(' ')
RT90_found=0;
RT50_found=0;
RT90_ind3=0;
RT50_ind3=0;
for ii=1:length(Trace_3Y),
    [mm3 mm_index3] =max(Trace_3Y);
    if ii > mm_index3 && (Trace_3Y(ii)/mm3)< .50 && RT50_found == 0; 
        RT50_ind3=ii;
        RT50_found=1;
    end
    if ii > mm_index3 && (Trace_3Y(ii)/mm3)< .10 && RT90_found == 0; 
        RT90_ind3=ii;
        RT90_found=1;
    end
    
end
disp(['Time to 50% Relaxation: ' num2str(Trace_3X(RT50_ind3)) 's for a pCa of ' num2str(pCa_Range(3)) ' at ' num2str(Pulse_Width) 'ms Twitch'])
disp(['Time to 90% Relaxation: ' num2str(Trace_3X(RT90_ind3)) 's for a pCa of ' num2str(pCa_Range(3)) ' at ' num2str(Pulse_Width) 'ms Twitch'])
%%

%Normalized Force/Time
subplot(3,1,2)
hold on;
[maximum max_idx] = max(Trace_1Y);
NormMax=mean(Trace_1Y(max_idx-3:max_idx+5));
% [maximum2 max_idx2] = max(Trace_2Y);
% [maximum3 max_idx3] = max(Trace_3Y);
plot(Trace_1X, Trace_1Y/NormMax,'k-'),hold on
plot(Trace_2X, Trace_2Y/NormMax,'r-'),hold on
plot(Trace_3X, Trace_3Y/NormMax,'g-'),hold on
xlabel('Time(s)')
ylabel('Normalized Force(pN)')
%axis([-.1 1.5 0 1],'on')
axis([-.1 2.0 0 1.2],'on')
set(gca,'YTick',0:.2:1)
set(gca,'XTick',0:.1:2)

%Calcium Profile/Time
subplot(3,1,3)
hold on;
plot(Trace_1X,Ca_level{1},'k-'),hold on
plot(Trace_2X,Ca_level{2},'r-'),hold on
plot(Trace_3X,Ca_level{3},'g-'),hold on
xlabel('Time(s)')
ylabel('pCa')
%axis([-.1 1.5 3 7],'on')
axis([-.1 2.0 3 7],'on')
set(gca, 'ydir', 'reverse')

%% Ton/Events/XB Force dot diagrams
clf(figure(3))
subplot(2,2,1)
hold on;
errorbar(pCa_Range(1),nanmean(ton_ms{1}), nanstd(ton_ms{1}), 'ko', 'markerfacecolor', 'k'), hold on
errorbar(pCa_Range(2),nanmean(ton_ms{2}), nanstd(ton_ms{2}), 'ro', 'markerfacecolor', 'r'), hold on
errorbar(pCa_Range(3),nanmean(ton_ms{3}), nanstd(ton_ms{3}), 'go', 'markerfacecolor', 'g'), hold on
set(gca, 'box', 'off')
ylabel('t_{on} (ms)')
axis([3 6 0 25],'on')
set(gca, 'xdir', 'reverse')
Title('t_{on} All')

subplot(2,2,2)
hold on;
errorbar(pCa_Range(1),nanmean(ton_Full{1}), nanstd(ton_Full{1}),'kd', 'markerfacecolor', 'k'), hold on
errorbar(pCa_Range(2),nanmean(ton_Full{2}), nanstd(ton_Full{2}), 'rd', 'markerfacecolor', 'r'), hold on
errorbar(pCa_Range(3),nanmean(ton_Full{3}), nanstd(ton_Full{3}), 'gd', 'markerfacecolor', 'g'), hold on
set(gca, 'box', 'off')
ylabel('t_{on} (ms)')
axis([3 6 0 25],'on')
set(gca, 'xdir', 'reverse')
Title('t_{on} Full')

subplot(2,2,3)
errorbar(pCa_Range(1),nanmean(ton_Back{1}), nanstd(ton_Back{1}), 'kx', 'markerfacecolor', 'k'), hold on
errorbar(pCa_Range(2),nanmean(ton_Back{2}), nanstd(ton_Back{2}), 'rx', 'markerfacecolor', 'r'), hold on
errorbar(pCa_Range(3),nanmean(ton_Back{3}), nanstd(ton_Back{3}), 'gx', 'markerfacecolor', 'g'), hold on
set(gca, 'box', 'off')
ylabel('t_{on} (ms)')
axis([3 6 0 25],'on')
set(gca, 'xdir', 'reverse')
Title('t_{on} Back')
%%

% subplot(2,2,2)
% errorbar(1, mean(N_raw)/mean(N_raw), std(ton_raw)/mean(N_raw), 'ko', 'markerfacecolor', 'k'), hold on
% errorbar(2, mean(Pct_full), std(Pct_full), 'ko', 'markerfacecolor', 'g'), hold on
% errorbar(3, mean(Pct_back), std(Pct_back), 'ko', 'markerfacecolor', 'm');
% set(gca, 'box', 'off')
% ylabel('Num. Events (percent)')
% xlim([0,4])
% ylim([0,1.2])
% set(gca, 'xtick', [1,2,3])
% set(gca, 'xticklabel', {'All', 'Full', 'Back Exit'})
% title(['N Total/trace ~' num2str(mean(N_raw)) '     ' date])

%% XB Force for All,Full, and Back Exit Cycles for 3 different pCa values
clf(figure(4))
subplot(2,2,1)
hold on;
errorbar(pCa_Range(1),nanmean(XBF_Raw{1}), nanstd(XBF_Raw{1}), 'ko', 'markerfacecolor', 'k'), hold on
errorbar(pCa_Range(2),nanmean(XBF_Raw{2}), nanstd(XBF_Raw{2}), 'ro', 'markerfacecolor', 'r');
errorbar(pCa_Range(3),nanmean(XBF_Raw{3}), nanstd(XBF_Raw{3}), 'go', 'markerfacecolor', 'g');
set(gca, 'box', 'off')
ylabel('XB Force (pN) Raw')
xlim([3,7])
ylim([-10,5])
title('Force from All Cycle XBs')
set(gca, 'xdir', 'reverse')
set(gca, 'xtick', pCa_Range)

subplot(2,2,2)
hold on;
errorbar(pCa_Range(1),nanmean(XBF_Full{1}), nanstd(XBF_Full{1}), 'kd', 'markerfacecolor', 'k'), hold on
errorbar(pCa_Range(2),nanmean(XBF_Full{2}), nanstd(XBF_Full{2}), 'rd', 'markerfacecolor', 'r');
errorbar(pCa_Range(3),nanmean(XBF_Full{3}), nanstd(XBF_Full{3}), 'gd', 'markerfacecolor', 'g');
set(gca, 'box', 'off')
ylabel('XB Force (pN) Full')
xlim([3,7])
ylim([-10,5])
title('Force from Full Cycle XBs')
set(gca, 'xdir', 'reverse')
set(gca, 'xtick', pCa_Range)

subplot(2,2,3)
hold on;
errorbar(pCa_Range(1),nanmean(XBF_Back{1}), nanstd(XBF_Back{1}), 'kx', 'markerfacecolor', 'k'), hold on
errorbar(pCa_Range(2),nanmean(XBF_Back{2}), nanstd(XBF_Back{2}), 'rx', 'markerfacecolor', 'r');
errorbar(pCa_Range(3),nanmean(XBF_Back{3}), nanstd(XBF_Back{3}), 'gx', 'markerfacecolor', 'g');
set(gca, 'box', 'off')
ylabel('XB Force (pN) Back Exit')
title('Force from Back Exit Cycle XBs')
xlim([3,7])
ylim([-10,5])
set(gca, 'xdir', 'reverse')
set(gca, 'xtick', pCa_Range)

%% XBF ENERGY for All,Full, and Back Exit Cycles for 3 different pCa values
clf(figure(5))
subplot(2,2,1)
hold on;
errorbar(pCa_Range(1),nanmean(XBF_EnRaw{1}), nanstd(XBF_EnRaw{1}), 'ko', 'markerfacecolor', 'k'), hold on
errorbar(pCa_Range(2),nanmean(XBF_EnRaw{2}), nanstd(XBF_EnRaw{2}), 'ro', 'markerfacecolor', 'r');
errorbar(pCa_Range(3),nanmean(XBF_EnRaw{3}), nanstd(XBF_EnRaw{3}), 'go', 'markerfacecolor', 'g');
set(gca, 'box', 'off')
ylabel('XB Energy (pN nm) Raw')
xlim([3,7])
ylim([0,100])
title('Energy from All Cycle XBs')
set(gca, 'xdir', 'reverse')
set(gca, 'xtick', pCa_Range)

subplot(2,2,2)
hold on;
errorbar(pCa_Range(1),nanmean(XBF_EnFull{1}), nanstd(XBF_EnFull{1}), 'kd', 'markerfacecolor', 'k'), hold on
errorbar(pCa_Range(2),nanmean(XBF_EnFull{2}), nanstd(XBF_EnFull{2}), 'rd', 'markerfacecolor', 'r');
errorbar(pCa_Range(3),nanmean(XBF_EnFull{3}), nanstd(XBF_EnFull{3}), 'gd', 'markerfacecolor', 'g');
set(gca, 'box', 'off')
ylabel('XB Energy (pN nm) Full')
xlim([3,7])
ylim([0,100])
Title('Energy from Full Cycle XBs')
set(gca, 'xdir', 'reverse')
set(gca, 'xtick', pCa_Range)

subplot(2,2,3)
hold on;
errorbar(pCa_Range(1),nanmean(XBF_EnBack{1}), nanstd(XBF_EnBack{1}), 'kx', 'markerfacecolor', 'k'), hold on
errorbar(pCa_Range(2),nanmean(XBF_EnBack{2}), nanstd(XBF_EnBack{2}), 'rx', 'markerfacecolor', 'r');
errorbar(pCa_Range(3),nanmean(XBF_EnBack{3}), nanstd(XBF_EnBack{3}), 'gx', 'markerfacecolor', 'g');
set(gca, 'box', 'off')
ylabel('XB Energy (pN nm) Back Exit')
xlim([3,7])
ylim([0,10])
Title('Energy from Back Exit XBs')
set(gca, 'xdir', 'reverse')
set(gca, 'xtick', pCa_Range)
%%
subplot(2,2,4)
errorbar(pCa_Range(1),nanmean(XBF_Enraw{1}), nanstd(XBF_raw{1}), 'ko', 'markerfacecolor', 'k'), hold on
errorbar(pCa_Range(2),nanmean(XBF_Enraw{2}), nanstd(XBF_raw{2}), 'ro', 'markerfacecolor', 'r');
errorbar(pCa_Range(3),nanmean(XBF_raw{3}), nanstd(XBF_raw{3}), 'go', 'markerfacecolor', 'g');
set(gca, 'box', 'off')
ylabel('XB Force (pN)')
xlim([3,7])
ylim([-10,5])
set(gca, 'xdir', 'reverse')
set(gca, 'xtick', pCa_Range)
errorbar(1, mean(XBNRG_raw), std(XBNRG_raw), 'ko', 'markerfacecolor', 'k');, hold on
errorbar(2, mean(XBNRG_full), std(XBNRG_full), 'ko', 'markerfacecolor', 'g');
errorbar(3, mean(XBNRG_back), std(XBNRG_back), 'ko', 'markerfacecolor', 'm');
set(gca, 'box', 'off')
ylabel('XB Energy (pN nm)')
xlim([0,4])
ylim([0,100])
set(gca, 'xtick', [1,2,3])
set(gca, 'xticklabel', {'All', 'Full', 'Back Exit'})


