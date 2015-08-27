function plotSS_ROIs

pCa =  num2str(4.0,'%3.2f');
RateRange = [0.1 0.25 0.5 1.0]; %-

SL = int2str(1210);
ShortenLength = int2str(-100);

ROIRange = ['A' 'B' 'C' 'D' 'E'];
%ROIRange = ['A' 'C' 'E'];
Muscle_Type_Range = {'N2B'};
colors = ['k'; 'b'; 'r'; 'g'; 'm'];
plottype = '.';

%clf(figure(1))
figure('Name','pCa = 4.00');


plot1 = []; plot2 = []; plot3 = [];
plot4 = []; plot5 = []; plot6 = [];
plot7 = []; plot8 = [];

for k = 1:length(Muscle_Type_Range)
    Muscle_Type = Muscle_Type_Range{k};
    for j = 1:length(ROIRange)
        ROI = ROIRange(j);
        
        subplot(3,3,2)
        xlim([0 1.1])
        hold on;
        xlabel('Rate')
        ylabel('Force(pN)')
        
        Values = zeros(1,length(RateRange));
        STDev = zeros(1,length(RateRange));
        for i = 1:length(RateRange)
            Rate = num2str(RateRange(i),'%3.2f');
            
            foldername = ['LengthChanges_Jan_Feb_2015/Muscle=', Muscle_Type, '_Rate=-',Rate,'_HalfSl=',SL,'/'];
            filepath = ['DataFiles/',foldername,'SSData_pCa_',pCa,'_ROI_',ROI,'.txt'];
            
            SSData = importdata(filepath);
            Values(i) = mean(SSData.data(:,3));
            STDev(i) = std(SSData.data(:,3));
            
            %plot(RateRange(i),mean(SSData.data(:,3)), [colors(j) plottype], 'MarkerSize', 15);
            %errorbar(RateRange(i),mean(SSData.data(:,3))),std(SSData.data(:,3));
        end
        plot1 = [plot1 plot(RateRange,Values, [colors(j) plottype], 'MarkerSize', 15)];
        legend(plot1,'1210 Before', 'Shortening','1110', 'Lengthening','1210 After', 'Location', 'NorthEast');
        errorbar(RateRange,Values,STDev);
        title([Muscle_Type ', pCa = ' pCa ',  SL = ' SL ShortenLength]);
        
        
        
        subplot(3,3,3)
        xlim([0 1.1])
        hold on;
        xlabel('Rate')
        ylabel('FA')
        
        Values = zeros(1,length(RateRange));
        STDev = zeros(1,length(RateRange));
        for i = 1:length(RateRange)
            Rate = num2str(RateRange(i),'%3.2f');
            
            foldername = ['LengthChanges_Jan_Feb_2015/Muscle=', Muscle_Type, '_Rate=-',Rate,'_HalfSl=',SL,'/'];
            filepath = ['DataFiles/',foldername,'SSData_pCa_',pCa,'_ROI_',ROI,'.txt'];
            
            SSData = importdata(filepath);
            Values(i) = mean(SSData.data(:,15));
            STDev(i) = std(SSData.data(:,15));
            
            %plot(RateRange(i),mean(SSData.data(:,15)), [colors(j,1) colors(j,2)], 'MarkerSize', 15);
        end
        plot2 = [plot2 plot(RateRange,Values, [colors(j) plottype], 'MarkerSize', 15)];
        errorbar(RateRange,Values,STDev);
        title([Muscle_Type ', pCa = ' pCa ',  SL = ' SL ShortenLength]);
        %legend(plot2,'1510 Before', 'Shortening','1410', 'Lengthening','1510 After', 'Location', 'NorthEast');
        
        
        subplot(3,3,4)
        xlim([0 1.1])
        hold on;
        xlabel('Rate')
        ylabel('xB')
        
        Values = zeros(1,length(RateRange));
        STDev = zeros(1,length(RateRange));
        for i = 1:length(RateRange)
            Rate = num2str(RateRange(i),'%3.2f');
            
            foldername = ['LengthChanges_Jan_Feb_2015/Muscle=', Muscle_Type, '_Rate=-',Rate,'_HalfSl=',SL,'/'];
            filepath = ['DataFiles/',foldername,'SSData_pCa_',pCa,'_ROI_',ROI,'.txt'];
            
            SSData = importdata(filepath);
            Values(i) = mean(SSData.data(:,19));
            STDev(i) = std(SSData.data(:,19));
            
            %plot(RateRange(i),mean(SSData.data(:,19)), [colors(j,1) colors(j,2)], 'MarkerSize', 15);
        end
        plot3 = [plot3 plot(RateRange,Values, [colors(j) plottype], 'MarkerSize', 15)];
        errorbar(RateRange,Values,STDev);
        title([Muscle_Type ', pCa = ' pCa ',  SL = ' SL ShortenLength]);
        %legend(plot3, '1510 Before', 'Shortening','1410', 'Lengthening','1510 After', 'Location', 'NorthEast');
        
        
        subplot(3,3,5)
        xlim([0 1.1])
        hold on;
        xlabel('Rate')
        ylabel('Force/xB')
        
        Values = zeros(1,length(RateRange));
        STDev = zeros(1,length(RateRange));
        for i = 1:length(RateRange)
            Rate = num2str(RateRange(i),'%3.2f');
            
            foldername = ['LengthChanges_Jan_Feb_2015/Muscle=', Muscle_Type, '_Rate=-',Rate,'_HalfSl=',SL,'/'];
            filepath = ['DataFiles/',foldername,'SSData_pCa_',pCa,'_ROI_',ROI,'.txt'];
            
            SSData = importdata(filepath);
            Values(i) = mean(SSData.data(:,3))/(mean(SSData.data(:,19))*720); %720 = number cross bridges
            STDev(i) = std(mean(SSData.data(:,3))/mean(SSData.data(:,19)));
            
            %plot(RateRange(i),mean(SSData.data(:,3))/mean(SSData.data(:,19)), [colors(j,1) colors(j,2)], 'MarkerSize', 15);
        end
        plot4 = [plot4 plot(RateRange,Values, [colors(j) plottype], 'MarkerSize', 15)];
        errorbar(RateRange,Values,STDev);
        title([Muscle_Type ', pCa = ' pCa ',  SL = ' SL ShortenLength]);
        %legend(plot4, '1510 Before', 'Shortening','1410', 'Lengthening','1510 After', 'Location', 'NorthEast');
        
        subplot(3,3,6)
        xlim([0 1.1])
        hold on;
        xlabel('Rate')
        ylabel('Ton Full')
        
        Values = zeros(1,length(RateRange));
        STDev = zeros(1,length(RateRange));
        for i = 1:length(RateRange)
            Rate = num2str(RateRange(i),'%3.2f');
            
            foldername = ['LengthChanges_Jan_Feb_2015/Muscle=', Muscle_Type, '_Rate=-',Rate,'_HalfSl=',SL,'/'];
            filepath = ['DataFiles/',foldername,'SS_ton_AVG_pCa_',pCa,'_ROI_',ROI,'.txt'];
            
            SSData = importdata(filepath);
            Values(i) = SSData.data(:,11);
            STDev(i) = SSData.data(:,12);
            
            %plot(RateRange(i),mean(SSData.data(:,3))/mean(SSData.data(:,19)), [colors(j,1) colors(j,2)], 'MarkerSize', 15);
        end
        plot5 = [plot5 plot(RateRange,Values, [colors(j) plottype], 'MarkerSize', 15)];
        errorbar(RateRange,Values,STDev);
        title([Muscle_Type ', pCa = ' pCa ',  SL = ' SL ShortenLength]);
        %legend(plot5, '1510 Before', 'Shortening','1410', 'Lengthening','1510 After', 'Location', 'NorthEast');
        
        subplot(3,3,7)
        xlim([0 1.1])
        hold on;
        xlabel('Rate')
        ylabel('% Ton Full')
        
        Values = zeros(1,length(RateRange));
        STDev = zeros(1,length(RateRange));
        for i = 1:length(RateRange)
            Rate = num2str(RateRange(i),'%3.2f');
            
            foldername = ['LengthChanges_Jan_Feb_2015/Muscle=', Muscle_Type, '_Rate=-',Rate,'_HalfSl=',SL,'/'];
            filepath = ['DataFiles/',foldername,'SS_ton_AVG_pCa_',pCa,'_ROI_',ROI,'.txt'];
            
            SSData = importdata(filepath);
            Values(i) = SSData.data(:,13);
            STDev(i) = SSData.data(:,14);
            
            %plot(RateRange(i),mean(SSData.data(:,3))/mean(SSData.data(:,19)), [colors(j,1) colors(j,2)], 'MarkerSize', 15);
        end
        plot6 = [plot6 plot(RateRange,Values, [colors(j) plottype], 'MarkerSize', 15)];
        errorbar(RateRange,Values,STDev);
        title([Muscle_Type ', pCa = ' pCa ',  SL = ' SL ShortenLength]);
        %legend(plot6,'1510 Before', 'Shortening','1410', 'Lengthening','1510 After', 'Location', 'NorthEast');
        
        
        subplot(3,3,8)
        xlim([0 1.1])
        hold on;
        xlabel('Rate')
        ylabel('Ton Back')
        
        Values = zeros(1,length(RateRange));
        STDev = zeros(1,length(RateRange));
        for i = 1:length(RateRange)
            Rate = num2str(RateRange(i),'%3.2f');
            
            foldername = ['LengthChanges_Jan_Feb_2015/Muscle=', Muscle_Type, '_Rate=-',Rate,'_HalfSl=',SL,'/'];
            filepath = ['DataFiles/',foldername,'SS_ton_AVG_pCa_',pCa,'_ROI_',ROI,'.txt'];
            
            SSData = importdata(filepath);
            Values(i) = SSData.data(:,15);
            STDev(i) = SSData.data(:,16);
            
            %plot(RateRange(i),mean(SSData.data(:,3))/mean(SSData.data(:,19)), [colors(j,1) colors(j,2)], 'MarkerSize', 15);
        end
        plot7 = [plot7 plot(RateRange,Values, [colors(j) plottype], 'MarkerSize', 15)];
        errorbar(RateRange,Values,STDev);
        title([Muscle_Type ', pCa = ' pCa ',  SL = ' SL ShortenLength]);
        %legend(plot7,'1510 Before', 'Shortening','1410', 'Lengthening','1510 After', 'Location', 'NorthEast');
        
        
        subplot(3,3,9)
        xlim([0 1.1])
        hold on;
        xlabel('Rate')
        ylabel('Ton (msec)')
        
        Values = zeros(1,length(RateRange));
        STDev = zeros(1,length(RateRange));
        for i = 1:length(RateRange)
            Rate = num2str(RateRange(i),'%3.2f');
            
            foldername = ['LengthChanges_Jan_Feb_2015/Muscle=', Muscle_Type, '_Rate=-',Rate,'_HalfSl=',SL,'/'];
            filepath = ['DataFiles/',foldername,'SS_ton_AVG_pCa_',pCa,'_ROI_',ROI,'.txt'];
            
            SSData = importdata(filepath);
            Values(i) = SSData.data(:,7);
            STDev(i) = SSData.data(:,8);
            
            %plot(RateRange(i),mean(SSData.data(:,3))/mean(SSData.data(:,19)), [colors(j,1) colors(j,2)], 'MarkerSize', 15);
        end
        plot8 = [plot8 plot(RateRange,Values, [colors(j) plottype], 'MarkerSize', 15)];
        errorbar(RateRange,Values,STDev);
        title([Muscle_Type ', pCa = ' pCa ',  SL = ' SL ShortenLength]);
        %legend(plot8,'1510 Before', 'Shortening','1410', 'Lengthening','1510 After', 'Location', 'NorthEast');
        
    end
end
