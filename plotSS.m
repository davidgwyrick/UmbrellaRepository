function [RateMatrix] = plotSS(pCaRange)

RateMatrix = zeros(size(pCaRange),4);

for i = 1:length(pCaRange)
    pCa = num2str(pCaRange(i),'%3.2f');
    StaticName=strcat('SSData_pCa_',pCa,'.txt');
    DynamicName=strcat('SSData_pCa_',pCa,'_LengthChange.txt');
    StaticFile=importdata(StaticName);
    DynamicFile=importdata(DynamicName);
    RateMatrix(i,1) = pCaRange(i);
    RateMatrix(i,2) = mean(StaticFile.data(:,3));
    RateMatrix(i,3) = mean(DynamicFile.data(:,3));
    RateMatrix(i,4) = (mean(StaticFile.data(:,3)))-(mean(DynamicFile.data(:,3)));
end


% clf(figure(1))
% subplot(3,1,1)
% hold on;
% 
% xlabel('Time(s)')
% ylabel('Force(pN)')

% for i = 1:length(pCaRange)
%     
%     pCa = num2str(pCaRange(i),'%3.2f');
%     filename=strcat('SSData_pCa_',pCa,'.txt');
%     
%     TimeSeries=importdata(filename);
%     
%     plot(TimeSeries.data(:,1),TimeSeries.data(:,2), 'k-')
%     plot(TimeSeries.data(:,1),TimeSeries.data(:,3), 'r-')
%     
% end
