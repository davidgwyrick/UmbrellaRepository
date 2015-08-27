%Takes in a pCa number or [numbers] as an agrument, then opens the corresponding
%text files in the working directory and plots the Force/Time and FA/Time relationships.
%
%Example: plotTS([4 5.25 6.5]) will open and plot these files together:
%   TimeSeriesAvg_pCa_4.00.txt
%   TimeSeriesAvg_pCa_5.25.txt
%   TimeSeriesAvg_pCa_6.50.txt

function plotTSwyrick(pCaRange,Outdir)

%% David Plotting
TopDir = pwd;
cd(Outdir)
clf(figure(1))
subplot(4,1,1)
hold on;

xlabel('Time(s)')
ylabel('Force(pN)')

    
    pCa = num2str(pCaRange(1),'%3.2f');
    filename=strcat('TimeSeriesAvg_pCa_',pCa,'.txt');
    
    TimeSeries=importdata(filename);
    %Top Plot of 4
    plot(TimeSeries.data(:,1),TimeSeries.data(:,2), 'k-')
    plot(TimeSeries.data(:,1),TimeSeries.data(:,3), 'r-')
   
    %Plot 2 of 4
    subplot(4,1,2)
    hold on;

    xlabel('Time(s)')
    ylabel('pCa')
    set(gca, 'ydir', 'reverse')
    plot(TimeSeries.data(:,1),TimeSeries.data(:,12), 'k-')
    
    %Plot 3 of 4
    subplot(4,1,3)
    hold on;
    
    xlabel('Time(s)')
    ylabel('TF Frac') %TF fraction available for binding
    plot(TimeSeries.data(:,1),TimeSeries.data(:,10), 'r-')
    
    %Plot 4 of 4
    subplot(4,1,4)
    hold on;

    xlabel('Time(s)')
    ylabel('Force/CrossBridge')
    plot(TimeSeries.data(:,1),TimeSeries.data(:,3)./(TimeSeries.data(:,5)*720), 'k-')
cd(TopDir)
end
%%




% clf(figure(1))
% subplot(4,1,1)
% hold on;
% 
% xlabel('Time(s)')
% ylabel('Force(pN)')
% 
% for i = 1:length(pCaRange)
%     
%     pCa = num2str(pCaRange(i),'%3.2f');
%     filename=strcat('TimeSeriesAvg_pCa_',pCa,'.txt');
%     
%     TimeSeries=importdata(filename);
%     
%     plot(TimeSeries.data(:,1),TimeSeries.data(:,2), 'k-')
%     plot(TimeSeries.data(:,1),TimeSeries.data(:,3), 'r-')
%     
% end
% 
% 
% % % clf(figure(2))
% % subplot(4,1,2)
% % hold on;
% % 
% % xlabel('Time(s)')
% % ylabel('Handle Position (nm)')
% % 
% % for i = 1:length(pCaRange)
% %     
% %     pCa = num2str(pCaRange(i),'%3.2f');
% %     filename=strcat('TimeSeriesAvg_pCa_',pCa,'.txt');
% %     
% %     TimeSeries=importdata(filename); %TimeSeriesAvg_pCa_4.00.txt
% %     
% %     plot(TimeSeries.data(:,1),TimeSeries.data(:,11), 'k-')
% %     
% % end
% 
% % clf(figure(2))
% %%
% subplot(4,1,2)
% hold on;
% 
% xlabel('Time(s)')
% ylabel('pCa')
% set(gca, 'ydir', 'reverse')
% 
% for i = 1:length(pCaRange)
%     
%     pCa = num2str(pCaRange(i),'%3.2f');
%     filename=strcat('TimeSeriesAvg_pCa_',pCa,'.txt');
%     
%     TimeSeries=importdata(filename); %TimeSeriesAvg_pCa_4.00.txt
%     plot(TimeSeries.data(:,1),TimeSeries.data(:,12), 'k-')
%     
% end
% %%
% 
% % clf(figure(3))
% subplot(4,1,3)
% hold on;
% 
% xlabel('Time(s)')
% ylabel('TF Frac')%TF fraction available for binding
% for i = 1:length(pCaRange)
%     
%     pCa = num2str(pCaRange(i),'%3.2f');
%     filename=strcat('TimeSeriesAvg_pCa_',pCa,'.txt');
%     
%     TimeSeries=importdata(filename); %TimeSeriesAvg_pCa_4.00.txt
%     
%     plot(TimeSeries.data(:,1),TimeSeries.data(:,10), 'r-')
% end 
% % clf(figure(3))
% subplot(4,1,4)
% hold on;
% 
% xlabel('Time(s)')
% ylabel('Force/CrossBridge')
% 
% for i = 1:length(pCaRange)
%     
%     pCa = num2str(pCaRange(i),'%3.2f');
%     filename=strcat('TimeSeriesAvg_pCa_',pCa,'.txt');
%     
%     TimeSeries=importdata(filename); %TimeSeriesAvg_pCa_4.00.txt
%     
%     plot(TimeSeries.data(:,1),TimeSeries.data(:,3)./(TimeSeries.data(:,5)*720), 'k-')
%     
% end
% 





