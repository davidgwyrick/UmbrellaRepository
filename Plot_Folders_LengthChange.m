%function Plot_Folders
clear all
%% MOD HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_rate = -0.25; %-0.1, -0.25, -0.5, or -1.0
Muscle_Type = 'Soleus';
%Muscle_Type = 'N2B';
%Trace 1
T1_SXB2 = 0.5; %0.5 or 1.0
T1_SXB3 = 0.1; %0.1 or 1.0
%Trace 2
T2_SXB2 = 1.0; %0.5 or 1.0
T2_SXB3 = 1.0; %0.1 or 1.0
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

pCa_Range = [4.0];
RateRange = [-0.1, -0.25, -0.5, -1.0]; 
ROI_Range = {'A','B','C','D','E'};
SXB2_Range = [1];
SXB3_Range = [1];

eta_range = [0.6850 0.8];
f1_range = [0.2800 0.5];
y0_range = [20 35]; 
SlopeL_range = [-100 -60];%, -80, -120]; 
SlopeR_range = [20 20];%, 24, 16];
atan_max_range = [100 50];%, 80, 120];
GausPeak_range = [2000 2000];%, 2400, 1600];

prcnt_Full = zeros(2,length(RateRange), length(ROI_Range));  %prcnt_Full @ (SXB2 index, SXB3 index, Rate, ROIs)
prcnt_Full_std = zeros(2,length(RateRange), length(ROI_Range));
ton_Full = zeros(2,length(RateRange), length(ROI_Range));
ton_Full_std = zeros(2,length(RateRange), length(ROI_Range));
ton_Back = zeros(2,length(RateRange), length(ROI_Range));
ton_Back_std = zeros(2,length(RateRange), length(ROI_Range));
ton_ms = zeros(2,length(RateRange), length(ROI_Range));
ton_ms_std = zeros(2,length(RateRange), length(ROI_Range));
Force = zeros(2,length(RateRange), length(ROI_Range));
Force_std = zeros(2,length(RateRange), length(ROI_Range));
FA = zeros(2,length(RateRange), length(ROI_Range));
FA_std = zeros(2,length(RateRange), length(ROI_Range));
XB = zeros(2,length(RateRange), length(ROI_Range));
XB_std = zeros(2,length(RateRange), length(ROI_Range));
ForceXB = zeros(2,length(RateRange), length(ROI_Range));
ForceXB_std = zeros(2,length(RateRange), length(ROI_Range));
TS_Time = {{},{}};  % TimeSeries Time @ (SXB2 index, SXB3 index, pCa index)
TS_Force = {{},{}}; % TimeSeries Force @ (SXB2 index, SXB3 index, pCa index)
Header_pos = {{},{}};

trace = 1;
    eta = eta_range(trace);
    f1 = f1_range(trace);
    y0 = y0_range(trace);
    SlopeL = SlopeL_range(trace);
    SlopeR = SlopeR_range(trace);
    atan_max = atan_max_range(trace);
    GausPeak = GausPeak_range(trace);
SXB2 = num2str(SXB2_Range(trace),'%1.2f');
SXB3 = num2str(SXB3_Range(trace),'%1.2f');



for iRate = 1:length(RateRange)
    Rate = num2str(RateRange(iRate));
%     for iSXB2 = 1:length(SXB2_Range)
%         SXB2 = num2str(SXB2_Range(iSXB2),'%1.2f');
%         for iSXB3 = 1:length(SXB3_Range)
%             SXB3 = num2str(SXB3_Range(iSXB3),'%1.2f');
            foldername = ['SXB2v3_LengthChange' filesep 'Muscle=',Muscle_Type,'_Rate=',Rate,'_HalfSl=1210_SXB2=',SXB2,'_SXB3=',SXB3,'/'];
            for ipCa = 1:length(pCa_Range)
                pCa = num2str(pCa_Range(ipCa),'%1.2f');
                for iROI = 1:length(ROI_Range)
                    ROI = ROI_Range{iROI};
                    filepath = ['J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex\DataFiles/',foldername,'SS_ton_AVG_pCa_',pCa,'_ROI_',ROI,'.txt'];
                    Data = importdata(filepath);
                    prcnt_Full(trace,iRate,iROI) = Data.data(13);
                    prcnt_Full_std(trace,iRate,iROI) = Data.data(14);
                    ton_Full(trace,iRate,iROI) = Data.data(11);
                    ton_Full_std(trace,iRate,iROI) = Data.data(12);
                    ton_Back(trace,iRate,iROI) = Data.data(15);
                    ton_Back_std(trace,iRate,iROI) = Data.data(16);
                    ton_ms(trace,iRate,iROI) = Data.data(7);
                    ton_ms_std(trace,iRate,iROI) = Data.data(8);
                    
                    filepath = ['J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex\DataFiles/',foldername,'SSData_pCa_',pCa,'_ROI_',ROI,'.txt'];
                    Data = importdata(filepath);
                    Force(trace,iRate,iROI) = mean(Data.data(:,3));
                    Force_std(trace,iRate,iROI) = std(Data.data(:,3));
                    FA(trace,iRate,iROI) = mean(Data.data(:,15));
                    FA_std(trace,iRate,iROI) = std(Data.data(:,15));
                    XB(trace,iRate,iROI) = mean(Data.data(:,19));
                    XB_std(trace,iRate,iROI) = std(Data.data(:,19));
                    
                    ForceXB(trace,iRate,iROI) = Force(trace,iRate,iROI)/(XB(trace,iRate,iROI)*720);
                    ForceXB_std(trace,iRate,iROI) = std(Data.data(:,3)./(Data.data(:,19)*720));
                    
                    filepath = ['J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex\DataFiles/',foldername,'TimeSeriesAvg_pCa_',pCa,'.txt'];
                    Data = importdata(filepath);
                    TS_Time(trace,iRate) = {Data.data(:,1)};
                    TS_Force(trace,iRate) = {Data.data(:,2)};
                    Header_pos(trace,iRate) = {Data.data(:,11)};
                end
            end
%         end
%     end
end



trace = 2;
    eta = eta_range(trace);
    f1 = f1_range(trace);
    y0 = y0_range(trace);
    SlopeL = SlopeL_range(trace);
    SlopeR = SlopeR_range(trace);
    atan_max = atan_max_range(trace);
    GausPeak = GausPeak_range(trace);


for iRate = 1:length(RateRange)
    Rate = num2str(RateRange(iRate));
%     for iSXB2 = 1:length(SXB2_Range)
%         SXB2 = num2str(SXB2_Range(iSXB2),'%1.2f');
%         for iSXB3 = 1:length(SXB3_Range)
%             SXB3 = num2str(SXB3_Range(iSXB3),'%1.2f');
            foldername = ['Kinetic_Rates' filesep 'Rate=', num2str(Rate),'SXB2=', SXB2,'_SXB3=', SXB3, '_eta=',num2str(eta,'%1.3f'),'_f1=', num2str(f1,'%1.3f'), '_y0=',num2str(y0,'%d'), '_SL=',num2str(SlopeL,'%d'), '_SR=',num2str(SlopeR,'%d'), '_atan=',num2str(atan_max,'%d'), '_GP=',num2str(GausPeak,'%d') filesep];
            for ipCa = 1:length(pCa_Range)
                pCa = num2str(pCa_Range(ipCa),'%1.2f');
                for iROI = 1:length(ROI_Range)
                    ROI = ROI_Range{iROI};
                    filepath = ['J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex\DataFiles/',foldername,'SS_ton_AVG_pCa_',pCa,'_ROI_',ROI,'.txt'];
                    Data = importdata(filepath);
                    prcnt_Full(trace,iRate,iROI) = Data.data(13);
                    prcnt_Full_std(trace,iRate,iROI) = Data.data(14);
                    ton_Full(trace,iRate,iROI) = Data.data(11);
                    ton_Full_std(trace,iRate,iROI) = Data.data(12);
                    ton_Back(trace,iRate,iROI) = Data.data(15);
                    ton_Back_std(trace,iRate,iROI) = Data.data(16);
                    ton_ms(trace,iRate,iROI) = Data.data(7);
                    ton_ms_std(trace,iRate,iROI) = Data.data(8);
                    
                    filepath = ['J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex\DataFiles/',foldername,'SSData_pCa_',pCa,'_ROI_',ROI,'.txt'];
                    Data = importdata(filepath);
                    Force(trace,iRate,iROI) = mean(Data.data(:,3));
                    Force_std(trace,iRate,iROI) = std(Data.data(:,3));
                    FA(trace,iRate,iROI) = mean(Data.data(:,15));
                    FA_std(trace,iRate,iROI) = std(Data.data(:,15));
                    XB(trace,iRate,iROI) = mean(Data.data(:,19));
                    XB_std(trace,iRate,iROI) = std(Data.data(:,19));
                    
                    ForceXB(trace,iRate,iROI) = Force(trace,iRate,iROI)/(XB(trace,iRate,iROI)*720);
                    ForceXB_std(trace,iRate,iROI) = std(Data.data(:,3)./(Data.data(:,19)*720));
                    
                    filepath = ['J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex\DataFiles/',foldername,'TimeSeriesAvg_pCa_',pCa,'.txt'];
                    Data = importdata(filepath);
                    TS_Time(trace,iRate) = {Data.data(:,1)};
                    TS_Force(trace,iRate) = {Data.data(:,2)};
                    Header_pos(trace,iRate) = {Data.data(:,11)};
                end
            end
%         end
%     end
end







%% Figure 1

Rate_index = find(RateRange == plot_rate); % index of desired pCa from pCa_Range (ie 2 = the second pCa value in the range). Optional: Create a pCa_index2 for comparison plots

% iT1_SXB2 = find(SXB2_Range == T1_SXB2);
% iT1_SXB3 = find(SXB3_Range == T1_SXB3);
% iT2_SXB2 = find(SXB2_Range == T2_SXB2);
% iT2_SXB3 = find(SXB3_Range == T2_SXB3);
iT1 = 1;
iT2 = 2;

Trace_1X = TS_Time{iT1,Rate_index};
Trace_1Y = TS_Force{iT1,Rate_index};
% Trace_1Label = ['SXB2/3 = ', num2str(SXB2_Range(iT1_SXB2),'%1.1f'), '/', num2str(SXB3_Range(iT1_SXB3),'%1.1f')];
Trace_1Label = 'Original Values';
Trace_2X = TS_Time{iT2,Rate_index};
Trace_2Y = TS_Force{iT2,Rate_index};
% Trace_2Label = ['SXB2/3 = ', num2str(SXB2_Range(iT2_SXB2),'%1.1f'), '/', num2str(SXB3_Range(iT2_SXB3),'%1.1f')];
Trace_2Label = 'New Values';

clf(figure(1))
subplot(3,1,1)
plot(Trace_1X, Trace_1Y, ...
    Trace_2X, Trace_2Y)
legend(Trace_1Label, ...
    Trace_2Label, ...
    'location', 'Northeast')
title(['TimeSeries pCa ', num2str(RateRange(Rate_index),'%1.2f'), ' Rate'])
xlabel('Time(s)')
ylabel('Force(pN)')

subplot(3,1,2)
% [maximum max_idx] = max(Trace_1Y);
% [maximum2 max_idx2] = max(Trace_2Y);
plot(Trace_1X, Trace_1Y/mean(Trace_1Y(1000:2000)), ...
    Trace_2X, Trace_2Y/mean(Trace_2Y(1000:2000)))
legend(Trace_1Label, ...
   Trace_2Label, ...
    'location', 'Northeast')
title(['Normalized TimeSeries pCa ', num2str(RateRange(Rate_index),'%1.2f'), ' Rate'])
xlabel('Time(s)')
ylabel('Force(pN)')

subplot(3,1,3)
plot(Trace_1X, Header_pos{trace,Rate_index})
title(['Rate = ', num2str(RateRange(Rate_index))])
xlabel('Time(s)')
ylabel('Header pos')

%% Figure 2
% SXB2 = 0.5;
% SXB3 = 0.1;

% iSXB2 = find(SXB2_Range == T1_SXB2);
% iSXB3 = find(SXB3_Range == T1_SXB3);

itrace = 1;

Col = {'k'; 'b'; 'r'; 'g'; 'm'};
clf(figure(2))

h2 = subplot(subplot(3,3,1));
set(gca,'visible','off')
text(0.3,1, sprintf('eta = %4.2f\nf1 = %4.2f\ny0 = %4.2f\nslopeL = %4.2f\natan MAX = %4.2f', eta_range(itrace), f1_range(itrace),y0_range(itrace),SlopeL_range(itrace),atan_max_range(itrace)), 'parent', h2);


subplot(3,3,2)
for n = 1:length(ROI_Range)
   errorbar(RateRange,squeeze(Force(itrace,:,n)),squeeze(Force_std(itrace,:,n)),['-',Col{n},'.'],'MarkerSize',15);
   hold on
end
legend('1210 Before', 'Shortening','1110', 'Lengthening','1210 After', 'Location', 'NorthEast');
xlim([-1.1 0])
ylim([0 2000])
xlabel('Rate')
ylabel('Force(pN)')
set(gca, 'xdir', 'reverse')
title([Muscle_Type ', pCa = ' pCa]);


subplot(3,3,3)
for n = 1:length(ROI_Range)
   errorbar(RateRange,squeeze(FA(itrace,:,n)),squeeze(FA_std(itrace,:,n)),['-',Col{n},'.'],'MarkerSize',15);
   hold on
end
xlim([-1.1 0])
ylim([0.990 1])
xlabel('Rate')
ylabel('FA')
set(gca, 'xdir', 'reverse')
title([Muscle_Type ', pCa = ' pCa]);

subplot(3,3,4)
for n = 1:length(ROI_Range)
   errorbar(RateRange,squeeze(XB(itrace,:,n)),squeeze(XB_std(itrace,:,n)),['-',Col{n},'.'],'MarkerSize',15);
   hold on
end
xlim([-1.1 0])
ylim([0.09 0.13])
xlabel('Rate')
ylabel('xB')
set(gca, 'xdir', 'reverse')
title([Muscle_Type ', pCa = ' pCa]);


subplot(3,3,5)
for n = 1:length(ROI_Range)
   errorbar(RateRange,squeeze(ForceXB(itrace,:,n)),squeeze(ForceXB_std(itrace,:,n)),['-',Col{n},'.'],'MarkerSize',15);
   hold on
end
xlim([-1.1 0])
ylim([0 25])
xlabel('Rate')
ylabel('Force/XB')
set(gca, 'xdir', 'reverse')
title([Muscle_Type ', pCa = ' pCa]);


subplot(3,3,6)
for n = 1:length(ROI_Range)
   errorbar(RateRange,squeeze(ton_Full(itrace,:,n)),squeeze(ton_Full_std(itrace,:,n)),['-',Col{n},'.'],'MarkerSize',15);
   hold on
end
xlim([-1.1 0])
ylim([0 25])
xlabel('Rate')
ylabel('Ton Full')
set(gca, 'xdir', 'reverse')
title([Muscle_Type ', pCa = ' pCa]);


subplot(3,3,7)
for n = 1:length(ROI_Range)
   errorbar(RateRange,squeeze(prcnt_Full(itrace,:,n)),squeeze(prcnt_Full_std(itrace,:,n)),['-',Col{n},'.'],'MarkerSize',15);
   hold on
end
xlim([-1.1 0])
ylim([0 0.4])
xlabel('Rate')
ylabel('% Ton Full')
set(gca, 'xdir', 'reverse')
title([Muscle_Type ', pCa = ' pCa]);


subplot(3,3,8)
for n = 1:length(ROI_Range)
   errorbar(RateRange,squeeze(ton_Back(itrace,:,n)),squeeze(ton_Back_std(itrace,:,n)),['-',Col{n},'.'],'MarkerSize',15);
   hold on
end
xlim([-1.1 0])
ylim([0 10])
xlabel('Rate')
ylabel('Ton Back')
set(gca, 'xdir', 'reverse')
title([Muscle_Type ', pCa = ' pCa]);


subplot(3,3,9)
for n = 1:length(ROI_Range)
   errorbar(RateRange,squeeze(ton_ms(itrace,:,n)),squeeze(ton_ms_std(itrace,:,n)),['-',Col{n},'.'],'MarkerSize',15);
   hold on
end
xlim([-1.1 0])
ylim([0 10])
xlabel('Rate')
ylabel('Ton (msec)')
set(gca, 'xdir', 'reverse')
title([Muscle_Type ', pCa = ' pCa]);

%% Figure 3
% SXB2 = 1.0;
% SXB3 = 1.0;

% iSXB2 = find(SXB2_Range == T2_SXB2);
% iSXB3 = find(SXB3_Range == T2_SXB3);

itrace = 2;

Col = {'k'; 'b'; 'r'; 'g'; 'm'};
clf(figure(3))

h2 = subplot(subplot(3,3,1));
set(gca,'visible','off')
text(0.3,1, sprintf('eta = %4.2f\nf1 = %4.2f\ny0 = %4.2f\nslopeL = %4.2f\natan MAX = %4.2f', eta_range(itrace), f1_range(itrace),y0_range(itrace),SlopeL_range(itrace),atan_max_range(itrace)), 'parent', h2);


subplot(3,3,2)
for n = 1:length(ROI_Range)
   errorbar(RateRange,squeeze(Force(itrace,:,n)),squeeze(Force_std(itrace,:,n)),['-',Col{n},'.'],'MarkerSize',15);
   hold on
end
legend('1210 Before', 'Shortening','1110', 'Lengthening','1210 After', 'Location', 'NorthEast');
xlim([-1.1 0])
ylim([0 2000])
xlabel('Rate')
ylabel('Force(pN)')
set(gca, 'xdir', 'reverse')
title([Muscle_Type ', pCa = ' pCa]);


subplot(3,3,3)
for n = 1:length(ROI_Range)
   errorbar(RateRange,squeeze(FA(itrace,:,n)),squeeze(FA_std(itrace,:,n)),['-',Col{n},'.'],'MarkerSize',15);
   hold on
end
xlim([-1.1 0])
ylim([0.990 1])
xlabel('Rate')
ylabel('FA')
set(gca, 'xdir', 'reverse')
title([Muscle_Type ', pCa = ' pCa]);

subplot(3,3,4)
for n = 1:length(ROI_Range)
   errorbar(RateRange,squeeze(XB(itrace,:,n)),squeeze(XB_std(itrace,:,n)),['-',Col{n},'.'],'MarkerSize',15);
   hold on
end
xlim([-1.1 0])
ylim([0.09 0.13])
xlabel('Rate')
ylabel('xB')
set(gca, 'xdir', 'reverse')
title([Muscle_Type ', pCa = ' pCa]);


subplot(3,3,5)
for n = 1:length(ROI_Range)
   errorbar(RateRange,squeeze(ForceXB(itrace,:,n)),squeeze(ForceXB_std(itrace,:,n)),['-',Col{n},'.'],'MarkerSize',15);
   hold on
end
xlim([-1.1 0])
ylim([0 25])
xlabel('Rate')
ylabel('Force/XB')
set(gca, 'xdir', 'reverse')
title([Muscle_Type ', pCa = ' pCa]);


subplot(3,3,6)
for n = 1:length(ROI_Range)
   errorbar(RateRange,squeeze(ton_Full(itrace,:,n)),squeeze(ton_Full_std(itrace,:,n)),['-',Col{n},'.'],'MarkerSize',15);
   hold on
end
xlim([-1.1 0])
ylim([0 25])
xlabel('Rate')
ylabel('Ton Full')
set(gca, 'xdir', 'reverse')
title([Muscle_Type ', pCa = ' pCa]);


subplot(3,3,7)
for n = 1:length(ROI_Range)
   errorbar(RateRange,squeeze(prcnt_Full(itrace,:,n)),squeeze(prcnt_Full_std(itrace,:,n)),['-',Col{n},'.'],'MarkerSize',15);
   hold on
end
xlim([-1.1 0])
ylim([0 0.4])
xlabel('Rate')
ylabel('% Ton Full')
set(gca, 'xdir', 'reverse')
title([Muscle_Type ', pCa = ' pCa]);


subplot(3,3,8)
for n = 1:length(ROI_Range)
   errorbar(RateRange,squeeze(ton_Back(itrace,:,n)),squeeze(ton_Back_std(itrace,:,n)),['-',Col{n},'.'],'MarkerSize',15);
   hold on
end
xlim([-1.1 0])
ylim([0 10])
xlabel('Rate')
ylabel('Ton Back')
set(gca, 'xdir', 'reverse')
title([Muscle_Type ', pCa = ' pCa]);


subplot(3,3,9)
for n = 1:length(ROI_Range)
   errorbar(RateRange,squeeze(ton_ms(itrace,:,n)),squeeze(ton_ms_std(itrace,:,n)),['-',Col{n},'.'],'MarkerSize',15);
   hold on
end
xlim([-1.1 0])
ylim([0 10])
xlabel('Rate')
ylabel('Ton (msec)')
set(gca, 'xdir', 'reverse')
title([Muscle_Type ', pCa = ' pCa]);