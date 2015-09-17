%% Ca Transient Profile
% 9/9/15
%
%% Load Land Data

filename=['HumanCai.txt'];
filepath=['J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\DataFiles\CaTransient\',filename];
humanCai=importdata(filepath);
filename=['mouseCai.txt'];
filepath=['J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\DataFiles\CaTransient\',filename];
mouseCai=importdata(filepath);
filename=['ratCai.txt'];
filepath=['J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\DataFiles\CaTransient\',filename];
ratCai=importdata(filepath);
filename=['fhuman.txt'];
filepath=['J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\DataFiles\CaTransient\',filename];
fhuman=importdata(filepath);
filename=['frat.txt'];
filepath=['J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\DataFiles\CaTransient\',filename];
frat=importdata(filepath);
filename=['fmouse.txt'];
filepath=['J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\DataFiles\CaTransient\',filename];
fmouse=importdata(filepath);

%% Load Wyrick Data
pulse_width_range=[100 160 300];
pCa_range=[4 5.6];
TS_Time=cell(1,9);
TS_Force=cell(1,9);
TS_Ca=cell(1,9);
Ca2state=cell(1,9);
index=1;
% for ii=1:length(pCa_range),
%     pCa=pCa_range(ii);
%     for jj=1:length(pulse_width_range),
%         pulse_width=pulse_width_range(jj);
%         foldername=[num2str(pulse_width),'msTwitch_pCapeak=',num2str(pCa)];
%         filepath=['J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\DataFiles\CaTransient\',foldername,'\TimeSeriesAvg_pCa_',num2str(pCa,'%1.2f'),'.txt'];
%         Data = importdata(filepath);
%         disp(foldername)
%         TS_Time(index) = {Data.data(:,1)};
%         TS_Force(index) = {Data.data(:,2)};
%         TS_Ca(index) = {Data.data(:,12)};
%         Ca2state(index) = {Data.data(:,10)};
%         index=index+1;
%     end
% end

filepath=['J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\DataFiles\CaTransient\Koff=100 RuOff=25 CaOff=100 KaTnI=100\Scaled_LandTwitchHuman_pCapeak=5.6\TimeSeriesAvg_pCa_5.60.txt'];
Data = importdata(filepath);
TS_Time(1)={Data.data(:,1)};
TS_Force(1) = {Data.data(:,2)};
TS_Ca(1) = {Data.data(:,12)};
Ca0state(1) = {Data.data(:,8)};
Ca1state(1) = {Data.data(:,9)};
Ca2state(1) = {Data.data(:,10)};

filepath=['J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\DataFiles\CaTransient\Koff=100 RuOff=25 CaOff=100 KaTnI=100\Scaled_LandTwitchMouse_pCapeak=5.6\TimeSeriesAvg_pCa_5.60.txt'];
Data = importdata(filepath);
TS_Time(2)={Data.data(:,1)};
TS_Force(2) = {Data.data(:,2)};
TS_Ca(2) = {Data.data(:,12)};
Ca0state(2) = {Data.data(:,8)};
Ca1state(2) = {Data.data(:,9)};
Ca2state(2) = {Data.data(:,10)};

filepath=['J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\DataFiles\CaTransient\Koff=100 RuOff=25 CaOff=100 KaTnI=100\Scaled_LandTwitchRat_pCapeak=5.6\TimeSeriesAvg_pCa_5.60.txt'];
Data = importdata(filepath);
TS_Time(3)={Data.data(:,1)};
TS_Force(3) = {Data.data(:,2)};
TS_Ca(3) = {Data.data(:,12)};
Ca0state(3) = {Data.data(:,8)};
Ca1state(3) = {Data.data(:,9)};
Ca2state(3) = {Data.data(:,10)};

filepath=['J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\DataFiles\CaTransient\Koff=150 RuOff=100 CaOff=100 KaTnI=10\Scaled_LandTwitchHuman_pCapeak=5.6\TimeSeriesAvg_pCa_5.60.txt'];
Data = importdata(filepath);
TS_Time(4)={Data.data(:,1)};
TS_Force(4) = {Data.data(:,2)};
TS_Ca(4) = {Data.data(:,12)};
Ca0state(4) = {Data.data(:,8)};
Ca1state(4) = {Data.data(:,9)};
Ca2state(4) = {Data.data(:,10)};

filepath=['J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\DataFiles\CaTransient\Koff=150 RuOff=100 CaOff=100 KaTnI=10\Scaled_LandTwitchMouse_pCapeak=5.6\TimeSeriesAvg_pCa_5.60.txt'];
Data = importdata(filepath);
TS_Time(5)={Data.data(:,1)};
TS_Force(5) = {Data.data(:,2)};
TS_Ca(5) = {Data.data(:,12)};
Ca0state(5) = {Data.data(:,8)};
Ca1state(5) = {Data.data(:,9)};
Ca2state(5) = {Data.data(:,10)};

filepath=['J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\DataFiles\CaTransient\Koff=150 RuOff=100 CaOff=100 KaTnI=10\Scaled_LandTwitchRat_pCapeak=5.6\TimeSeriesAvg_pCa_5.60.txt'];
Data = importdata(filepath);
TS_Time(6)={Data.data(:,1)};
TS_Force(6) = {Data.data(:,2)};
TS_Ca(6) = {Data.data(:,12)};
Ca0state(6) = {Data.data(:,8)};
Ca1state(6) = {Data.data(:,9)};
Ca2state(6) = {Data.data(:,10)};

filepath=['J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\DataFiles\CaTransient\Koff=150 RuOff=100 CaOff=100 KaTnI=100\Scaled_LandTwitchHuman_pCapeak=5.6\TimeSeriesAvg_pCa_5.60.txt'];
Data = importdata(filepath);
TS_Time(7)={Data.data(:,1)};
TS_Force(7) = {Data.data(:,2)};
TS_Ca(7) = {Data.data(:,12)};
Ca0state(7) = {Data.data(:,8)};
Ca1state(7) = {Data.data(:,9)};
Ca2state(7) = {Data.data(:,10)};

filepath=['J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\DataFiles\CaTransient\Koff=150 RuOff=100 CaOff=100 KaTnI=100\Scaled_LandTwitchMouse_pCapeak=5.6\TimeSeriesAvg_pCa_5.60.txt'];
Data = importdata(filepath);
TS_Force(8) = {Data.data(:,2)};
TS_Ca(8) = {Data.data(:,12)};
Ca0state(8) = {Data.data(:,8)};
Ca1state(8) = {Data.data(:,9)};
Ca2state(8) = {Data.data(:,10)};

filepath=['J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\DataFiles\CaTransient\Koff=150 RuOff=100 CaOff=100 KaTnI=100\Scaled_LandTwitchRat_pCapeak=5.6\TimeSeriesAvg_pCa_5.60.txt'];
Data = importdata(filepath);
TS_Force(9) = {Data.data(:,2)};
TS_Ca(9) = {Data.data(:,12)};
Ca0state(9) = {Data.data(:,8)};
Ca1state(9) = {Data.data(:,9)};
Ca2state(9) = {Data.data(:,10)};
%% 
% NSTEPS=1001
% CaWyHuman = (0.25*10^-7)*ones(NSTEPS,1);
% Ca_0 = 0.14; % 0.1399*10^(-6)molar Resting Diastolic Calcium concentration 
% Ca_max =.5; %0.3431*10^(-6)molar is pCa of 6.46 Ca Transient Max concentration at Time to Peak (48.2ms)
%  tao=DetermineTaoParam(600,6.3);
%  disp('Human')
%  disp(tao)
%   tao=80;
% for indy=1:NSTEPS 
%     CaWyHuman(indy) = Ca_0 + (Ca_max - Ca_0)*((indy-1)/tao)*exp(1-(indy/tao));
% end
% CaWyHuman=CaWyHuman';
% %%
% NSTEPS=168;
% CaWyMouse=ones(NSTEPS,1);
% Ca_0 = 0.199;
% Ca_max=.4965;
%  tao=DetermineTaoParam(55,6.3);
%   disp('Mouse')
%  disp(tao)
% tao=25;
% for indy=1:NSTEPS 
%     CaWyMouse(indy) = Ca_0 + (Ca_max - Ca_0)*((indy-1)/tao)*exp(1-(indy/tao));
% end
% CaWyMouse=CaWyMouse';
% %%
% NSTEPS=168;
% CaWyRat=ones(NSTEPS,1);
% Ca_0 = 0.237;
% Ca_max=1.423;
%  tao=DetermineTaoParam(65,5.85);
%   disp('Rat')
%  disp(tao)
% tao=25;
% for indy=1:NSTEPS 
%     CaWyRat(indy) = Ca_0 + (Ca_max - Ca_0)*((indy-1)/tao)*exp(1-(indy/tao));
% end
% CaWyRat=CaWyRat';
% %%


%% Setting up the Force, TF3 state, and Ca profile vectors
%100ms Twitch w/ peak pCa=4
Trace1Y=TS_Force{1};
TF1state1=Ca0state{1};
TF2state1=Ca1state{1};
TF3state1=Ca2state{1};
CaLevel1=(10^6).*10.^(TS_Ca{1});
f5= ones(1,7)/7;%Average over 7 data points or 7ms
avgTrace1Y=filtfilt(f5,1,Trace1Y);
ym1=max(avgTrace1Y);

% plot(TS_Time{1},Trace1Y,'-b'),hold on
% plot(TS_Time{1},avgTrace1Y,'--m');
% legend('Raw Force Data', 'Moving Average Zero Phase Shift')
% xlim([0 1]);
% xlabel('Time');
% ylabel('Force');

%160ms Twitch w/ peak pCa=4
Trace2Y=TS_Force{2};
TF1state2=Ca0state{2};
TF2state2=Ca1state{2};
TF3state2=Ca2state{2};
CaLevel2=(10^6).*10.^(TS_Ca{2});
avgTrace2Y=filtfilt(f5,1,Trace2Y);
ym2=max(avgTrace2Y);

%300ms Twitch w/ peak pCa=4
Trace3Y=TS_Force{3};
TF1state3=Ca0state{3};
TF2state3=Ca1state{3};
TF3state3=Ca2state{3};
CaLevel3=(10^6).*10.^(TS_Ca{3});
avgTrace3Y=filtfilt(f5,1,Trace3Y);
ym3=max(avgTrace3Y);

%100ms Twitch w/ peak pCa=5.6
Trace4Y=TS_Force{4};
TF1state4=Ca0state{4};
TF2state4=Ca1state{4};
TF3state4=Ca2state{4};
CaLevel4=(10^6).*10.^(-TS_Ca{4});
avgTrace4Y=filtfilt(f5,1,Trace4Y);
ym4=max(avgTrace4Y);

%160ms Twitch w/ peak pCa=5.6
Trace5Y=TS_Force{5};
TF1state5=Ca0state{5};
TF2state5=Ca1state{5};
TF3state5=Ca2state{5};
CaLevel5=(10^6).*10.^(-TS_Ca{5});
avgTrace5Y=filtfilt(f5,1,Trace5Y);
ym5=max(avgTrace5Y);

%300ms Twitch w/ peak pCa=5.6
Trace6Y=TS_Force{6};
TF1state6=Ca0state{6};
TF2state6=Ca1state{6};
TF3state6=Ca2state{6};
CaLevel6=(10^6).*10.^(-TS_Ca{6});
avgTrace6Y=filtfilt(f5,1,Trace6Y);
ym6=max(avgTrace6Y);

%Land's Scaled Human Ca Transient Run through Tanner Model
Trace7Y=TS_Force{7};
TF1state7=Ca0state{7};
TF2state7=Ca1state{7};
TF3state7=Ca2state{7};
CaLevel7=(10^6).*10.^(-TS_Ca{7});
avgTrace7Y=filtfilt(f5,1,Trace7Y);
ym7=max(avgTrace7Y);

%Land's Scaled Mouse Ca Transient Run through Tanner Model
Trace8Y=TS_Force{8};
TF1state8=Ca0state{8};
TF2state8=Ca1state{8};
TF3state8=Ca2state{8};
CaLevel8=(10^6).*10.^(-TS_Ca{8});
avgTrace8Y=filtfilt(f5,1,Trace8Y);
ym8=max(avgTrace8Y);

%Land's Scaled Rat Ca Transient Run through Tanner Model
Trace9Y=TS_Force{9};
TF1state9=Ca0state{9};
TF2state9=Ca1state{9};
TF3state9=Ca2state{9};
CaLevel9=(10^6).*10.^(-TS_Ca{9});
avgTrace9Y=filtfilt(f5,1,Trace9Y);
ym9=max(avgTrace9Y);

timeTrace2=(10^3).*TS_Time{1}; %For 2000 step simulation
timeTrace6=(10^3).*TS_Time{4}; %For 6000 step simulation

%% Ca Profile Curves
figure(1); clf; hold on; 
% plot(0:1001,[mouseCai(1:end-1) mouseCai(1:end-1) mouseCai(1:end-1) mouseCai(1:end-1) mouseCai(1:end-1) mouseCai(1:end-1)],'-','Color',[0.8 0.2 0.2],'LineWidth',2);
% plot(0:1001,[ratCai(1:end-1) ratCai(1:end-1) ratCai(1:end-1) ratCai(1:end-1) ratCai(1:end-1) ratCai(1:end-1)],'-','Color',[0.2 0.2 0.8],'LineWidth',2);
% plot(0:1000,humanCai,'-','Color',[0.2 0.8 0.2],'LineWidth',2);

plot(timeTrace6,CaLevel8,'-','Color',[0.8 0.2 0.2],'LineWidth',1.5),hold on
plot(timeTrace6,CaLevel9,'-','Color',[0.2 0.2 0.8],'LineWidth',1.5),hold on
plot(timeTrace6,CaLevel7,'-','Color',[0.2 0.8 0.2],'LineWidth',1.5),hold on
plot(timeTrace6,CaLevel4,'--k','LineWidth',1.5),hold on
plot(timeTrace6,CaLevel5,'--m','LineWidth',1.5),hold on
plot(timeTrace6,CaLevel6,'--c','LineWidth',1.5),hold on




xlim([0 500]); 
xlabel('Time (ms)'); ylabel('[Ca^{2+}]_i \mu{}M');
legend('Scaled Mouse', 'Scaled Rat', 'Scaled Human','100ms-pCa=5.6','160ms-pCa=5.6','300ms-pCa=5.6');

%% Not so messy averaged Curves
% 

% clf(figure(11));% Not so messy averaged Curves
% plot(timeTrace2,avgTrace1Y/ym1,'-k','LineWidth',1.5),hold on
% plot(timeTrace2,avgTrace2Y/ym2,'-g','LineWidth',1.5),hold on
% plot(timeTrace2,avgTrace3Y/ym3,'-b','LineWidth',1.5),hold on
% plot(timeTrace6,avgTrace4Y/ym4,'-r','LineWidth',1.5),hold on
% plot(timeTrace6,avgTrace5Y/ym5,'-c','LineWidth',1.5),hold on
% plot(timeTrace6,avgTrace6Y/ym6,'-m','LineWidth',1.5),hold on
% axis([0 1000 0 1],'on');
% xlabel('Time (ms)')
% ylabel('Normalized Force');
% legend('100ms pCaPeak=4','160ms pCaPeak=4','300ms pCaPeak=4','100ms pCaPeak=5.6','160ms pCaPeak=5.6','300ms pCaPeak=5.6');
% title('Not-so-Messy Moving Averaged Force-Time curves');
% %Messy non-averaged curves
% clf(figure(12)); 
% plot(TS_Time{1},Trace1Y/ym1,'-k'),hold on
% plot(TS_Time{2},Trace2Y/ym2,'-g'),hold on
% plot(TS_Time{3},Trace3Y/ym3,'-b'),hold on
% plot(TS_Time{4},Trace4Y/ym4,'-r'),hold on
% plot(TS_Time{5},Trace5Y/ym5,'-c'),hold on
% plot(TS_Time{6},Trace6Y/ym6,'-m'),hold on
% axis([0 1 0 1],'on');
% xlabel('Time (ms)')
% ylabel('Normalized Force');
% legend('100ms pCaPeak=4','160ms pCaPeak=4','300ms pCaPeak=4','100ms pCaPeak=5.6','160ms pCaPeak=5.6','300ms pCaPeak=5.6');
% title('Messy Non-Averaged Force-Time curves');
%% Normalized Forces (Land)
%K=6.9565;
fmMax=max(fmouse(end-500:end));
fhMax=max(fhuman(end-1000:end));
frMax=max(frat(end-500:end));
% HumanTrace=TS_Force{7};
% HTmax=max(HumanTrace);
% MouseTrace=TS_Force{8};
% MTmax=max(MouseTrace);
% RatTrace=TS_Force{9};
% RTmax=max(RatTrace);

clf(figure(2));
subplot(2,2,1)
hold on;
plot(0:500,fmouse(end-500:end)/fmMax,'--g','LineWidth',2),hold on
plot(timeTrace6,avgTrace2Y/ym2,'-r','LineWidth',1.5),hold on
%plot(timeTrace6,filtfilt(f5,1,TF3state2),':r','LineWidth',1.5),hold on
%plot(timeTrace6,avgTrace5Y/ym5,'-b','LineWidth',1.5),hold on

plot(timeTrace6,avgTrace8Y/ym8,'-b','LineWidth',1.5),hold on
plot(timeTrace6,filtfilt(f5,1,TF3state8),':k','LineWidth',1.5)
xlim([0 500]); 
xlabel('Time (ms)');
ylabel('Normalized Force');
title('Mouse');
legend('Land Model','Koff=CaOff=100 RuOff=25 KaTnI=100','Koff=150 RuOff=100 CaOff=100 KaTnI=100','TF3 Fraction');

subplot(2,2,2)
hold on;
plot(0:500,frat(end-500:end)/frMax,'--g','LineWidth',2),hold on
plot(timeTrace6,avgTrace3Y/ym3,'-r','LineWidth',1.5)
%plot(timeTrace6,filtfilt(f5,1,TF3state3),':r','LineWidth',1.5),hold on
%plot(timeTrace6,avgTrace6Y/ym6,'-b','LineWidth',1.5)

plot(timeTrace6,avgTrace9Y/ym9,'-b','LineWidth',1.5),hold on
plot(timeTrace6,filtfilt(f5,1,TF3state9),':k','LineWidth',1.5)
xlim([0 500]); 
xlabel('Time (ms)');
ylabel('Normalized Force');
title('Rat');
legend('Land Model','Koff=CaOff=100 RuOff=25 KaTnI=100','Koff=150 RuOff=100 CaOff=100 KaTnI=100','TF3 Fraction');

subplot(2,2,3)
hold on;
plot(0:1000,fhuman(end-1000:end)/fhMax,'--g','LineWidth',2),hold on
plot(timeTrace6,avgTrace1Y/ym1,'-r','LineWidth',1.5),hold on
%plot(timeTrace6,filtfilt(f5,1,TF3state1),':r','LineWidth',1.5),hold on
%plot(timeTrace6,avgTrace4Y/ym4,'-b','LineWidth',1.5),hold on

plot(timeTrace6,avgTrace7Y/ym7,'-b','LineWidth',1.5),hold on
plot(timeTrace6,filtfilt(f5,1,TF3state7),':k','LineWidth',1.5)
xlim([0 500]); 
xlabel('Time (ms)');
ylabel('Normalized Force');
title('Human');
legend('Land Model','Koff=CaOff=100 RuOff=25 KaTnI=100','Koff=150 RuOff=100 CaOff=100 KaTnI=100','TF3 Fraction');

%% Comparing TF states
clf(figure(13));
subplot(2,2,1)
hold on;

plot(timeTrace6,filtfilt(f5,1,TF3state2),'--r','LineWidth',1.5),hold on
plot(timeTrace6,filtfilt(f5,1,TF3state5),'--g','LineWidth',1.5)
plot(timeTrace6,filtfilt(f5,1,TF3state8),'--b','LineWidth',1.5)
xlim([0 500]); 
xlabel('Time (ms)');
ylabel('TF3 Fraction');
title('Mouse');
legend('Koff=CaOff=100 RuOff=25 KaTnI=100','Koff= 150 RuOff=100 CaOff=100 KaTnI=10','Koff=150 RuOff=100 CaOff=100 KaTnI=100');

subplot(2,2,2)
hold on;
plot(timeTrace6,filtfilt(f5,1,TF3state3),'--r','LineWidth',1.5),hold on
plot(timeTrace6,filtfilt(f5,1,TF3state6),'--g','LineWidth',1.5)
plot(timeTrace6,filtfilt(f5,1,TF3state9),'--b','LineWidth',1.5)
xlim([0 500]); 
xlabel('Time (ms)');
ylabel('TF3 Fraction');
title('Rat');
legend('Koff=CaOff=100 RuOff=25 KaTnI=100','Koff= 150 RuOff=100 CaOff=100 KaTnI=10','Koff=150 RuOff=100 CaOff=100 KaTnI=100');

subplot(2,2,3)
hold on;
plot(timeTrace6,filtfilt(f5,1,TF3state1),'--r','LineWidth',1.5),hold on
plot(timeTrace6,filtfilt(f5,1,TF3state4),'--g','LineWidth',1.5)
plot(timeTrace6,filtfilt(f5,1,TF3state7),'--b','LineWidth',1.5)
xlim([0 500]); 
xlabel('Time (ms)');
ylabel('TF3 Fraction');
title('Human');
legend('Koff=CaOff=100 RuOff=25 KaTnI=100','Koff= 150 RuOff=100 CaOff=100 KaTnI=10','Koff=150 RuOff=100 CaOff=100 KaTnI=100');
%%
plot(timeTrace2,avgTrace1Y/ym1,'-k','LineWidth',1.5),hold on
plot(timeTrace2,avgTrace2Y/ym2,'-g','LineWidth',1.5),hold on
plot(timeTrace2,avgTrace3Y/ym3,'-b','LineWidth',1.5),hold on
plot(timeTrace6,avgTrace4Y/ym4,'-r','LineWidth',1.5),hold on
plot(timeTrace6,avgTrace5Y/ym5,'-c','LineWidth',1.5),hold on
plot(timeTrace6,avgTrace6Y/ym6,'-m','LineWidth',1.5),hold on
xlim([0 500]); 
xlabel('Time (ms)');
ylabel('Normalized Force');
legend('Mouse', 'Rat', 'Human','100ms pCa=4','160ms pCa=4','300ms pCa=4');
  