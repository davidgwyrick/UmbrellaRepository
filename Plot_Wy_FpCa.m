%pCa_Range = [4.0, 4.5, 5.0, 5.5, 5.7, 5.8, 5.9, 6.0, 6.1, 6.25, 6.5, 7.0, 8.0];
pCa_Range= [4 5 5.5 5.75 6 6.5 7 8];
%RateRange = [0]; %-


prcnt_Full = zeros(1,length(pCa_Range));
ton_Full = zeros(1,length(pCa_Range));
ton_Back = zeros(1,length(pCa_Range));
ton_ms = zeros(1,length(pCa_Range));
Force = zeros(1,length(pCa_Range));
Force_std = zeros(1,length(pCa_Range));


%addpath('J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex')



foldername = ['SSdataTest\'];
for ipCa = 1:length(pCa_Range)
    pCa = num2str(pCa_Range(ipCa),'%1.2f');
    %filepath = ['J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex\DataFiles/',foldername,'SS_ton_AVG_pCa_',pCa,'_ROI_A.txt'];
    filepath = ['J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\DataFiles\',foldername,'SS_ton_AVG_pCa_',pCa,'_ROI_A.txt'];
    Data = importdata(filepath);
    prcnt_Full(1,ipCa) = Data.data(13);
    ton_Full(1,ipCa) = Data.data(11);
    ton_Back(1,ipCa) = Data.data(15);
    ton_ms(1,ipCa) = Data.data(7);

    filepath = ['J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\DataFiles\',foldername,'SSData_pCa_',pCa,'_ROI_A.txt'];
    Data = importdata(filepath);
    Force(1,ipCa) = mean(Data.data(:,3));
    Force_std(1,ipCa) = std(Data.data(:,3));
end
 
clf(figure(3))
plot(pCa_Range,Force)
set(gca,'XDir','Reverse')
xlabel('pCa');
ylabel('Force (pN)');

% clf(figure(3))
% plot(pCa_Range, squeeze(Force(1,1,:)), pCa_Range, squeeze(Force(2,2,:)))
% % plot(pCa_Range, squeeze(Force(1,3,:))), hold on
% % plot(pCa_Range, squeeze(Force(3,1,:))), hold on
% set(gca,'XDir','Reverse')
% xlabel('pCa');
% ylabel('Force (pN)');
% %legend('SXB2/3 = 0.5/0.1','SXB2/3 = 1/1',  'location', 'Northwest')

% 
% Force(:,:,5) = 0; 
% Force_std(:,:,5) = 0;
% pCa_Range(5) = 8.0;

[head1, coeff1]=Process_3ParamHill_Anal_v2(pCa_Range', squeeze(Force(1,1,:)));
[head2, coeff2]=Process_3ParamHill_Anal_v2(pCa_Range', squeeze(Force(2,2,:)));
% [head3, coeff3]=Process_3ParamHill_Anal_v2(pCa_Range', squeeze(Force(1,3,:)));
% [head4, coeff4]=Process_3ParamHill_Anal_v2(pCa_Range', squeeze(Force(3,1,:)));

fprintf('\nFor SXB2 of 0.5 and SXB3 of 0.1\nFmax = %4.2f\nHill Coeff = %4.2f\npCa50 = %4.2f [%4.2f to %4.2f]\n',coeff1(1,1), coeff1(2,1), coeff1(3,1), coeff1(3,4),coeff1(3,5))
fprintf('\nFor SXB2 of 1 and SXB3 of 1\nFmax = %4.2f\nHill Coeff = %4.2f\npCa50 = %4.2f [%4.2f to %4.2f]\n',coeff2(1,1), coeff2(2,1), coeff2(3,1), coeff2(3,4),coeff2(3,5))

text(50,50, sprintf('For SXB2 of 0.5 and SXB3 of 0.1\nFmax = %4.2f\nHill Coeff = %4.2f\npCa50 = %4.2f [%4.2f to %4.2f]\n\nFor SXB2 of 1 and SXB3 of 1\nFmax = %4.2f\nHill Coeff = %4.2f\npCa50 = %4.2f [%4.2f to %4.2f]\n',coeff1(1,1), coeff1(2,1), coeff1(3,1), coeff1(3,4),coeff1(3,5),coeff2(1,1), coeff2(2,1), coeff2(3,1), coeff2(3,4),coeff2(3,5)));

clf(figure(4))
subplot(2,2,1)
plot([4:0.02:9], (coeff1(1,1)./(1+10.^(coeff1(2,1)*([4:0.02:9] - coeff1(3,1))))),...
     [4:0.02:9], (coeff2(1,1)./(1+10.^(coeff2(2,1)*([4:0.02:9] - coeff2(3,1)))))), hold on
%     [4:0.02:9], (coeff3(1,1)./(1+10.^(coeff3(2,1)*([4:0.02:9] - coeff3(3,1))))),...
%     [4:0.02:9], (coeff4(1,1)./(1+10.^(coeff4(2,1)*([4:0.02:9] - coeff4(3,1)))))), hold on

legend('SXB2/3 = 0.5/0.1','SXB2/3 = 1/1',  'location', 'Northwest')%,'SXB2/3 = 0.1/0.5','SXB2/3 = 0.5/0.1',  'location', 'Northwest')
errorbar(pCa_Range, squeeze(Force(1,1,:)), squeeze(Force_std(1,1,:)),'bo')
errorbar(pCa_Range, squeeze(Force(2,2,:)), squeeze(Force_std(2,2,:)),'ko')
% errorbar(pCa_Range, squeeze(Force(1,3,:)), squeeze(Force_std(1,3,:)), 'ro')
% errorbar(pCa_Range, squeeze(Force(3,1,:)), squeeze(Force_std(3,1,:)), 'co')
title('Force/pCa')
ylabel('Force (pN)')
xlabel('pCa')
set(gca, 'xdir', 'reverse')
set(gca, 'box', 'off')

h2 = subplot(subplot(2,2,2));
set(gca,'visible','off')
text(0.2,0.4, sprintf('For SXB2 of 0.5 and SXB3 of 0.1\nFmax = %4.2f\nHill Coeff = %4.2f\npCa50 = %4.2f [%4.2f to %4.2f]\n\nFor SXB2 of 1 and SXB3 of 1\nFmax = %4.2f\nHill Coeff = %4.2f\npCa50 = %4.2f [%4.2f to %4.2f]\n',coeff1(1,1), coeff1(2,1), coeff1(3,1), coeff1(3,4),coeff1(3,5),coeff2(1,1), coeff2(2,1), coeff2(3,1), coeff2(3,4),coeff2(3,5)), 'parent', h2);


subplot(2,2,3)
plot([4:0.02:9], (coeff1(1,1)./(1+10.^(coeff1(2,1)*([4:0.02:9] - coeff1(3,1)))) / coeff2(1,1)), ...
     [4:0.02:9], (coeff2(1,1)./(1+10.^(coeff2(2,1)*([4:0.02:9] - coeff2(3,1)))) / coeff2(1,1))), hold on
%    [4:0.02:9], (coeff3(1,1)./(1+10.^(coeff3(2,1)*([4:0.02:9] - coeff3(3,1)))) / coeff1(1,1)), ...
%    [4:0.02:9], (coeff4(1,1)./(1+10.^(coeff4(2,1)*([4:0.02:9] - coeff4(3,1)))) / coeff1(1,1))), hold on

 
legend('SXB2/3 = 0.5/0.1','SXB2/3 = 1/1',  'location', 'Northwest')%,'SXB2/3 = 0.1/0.5','SXB2/3 = 0.5/0.1',  'location', 'Northwest')
errorbar(pCa_Range, (squeeze(Force(1,1,:)))/ Force(2,2,1), squeeze(Force_std(1,1,:))./Force(2,2,1), 'bo')
errorbar(pCa_Range, (squeeze(Force(2,2,:))/ Force(2,2,1)), squeeze(Force_std(2,2,:))./Force(2,2,1), 'ko')
% errorbar(pCa_Range, (squeeze(Force(1,3,:))/ Force(4,4,1)), squeeze(Force_std(1,3,:))./Force(4,4,1), 'ro')
% errorbar(pCa_Range, (squeeze(Force(3,1,:))/ Force(4,4,1)), squeeze(Force_std(3,1,:))./Force(4,4,1), 'co')
title('Relative Force/pCa')
ylabel('Relative Force (pN)')
xlabel('pCa')
set(gca, 'xdir', 'reverse')
set(gca, 'box', 'off')

subplot(2,2,4)
plot([4:0.02:9], (coeff1(1,1)./(1+10.^(coeff1(2,1)*([4:0.02:9] - coeff1(3,1)))) / coeff1(1,1)), ...
     [4:0.02:9], (coeff2(1,1)./(1+10.^(coeff2(2,1)*([4:0.02:9] - coeff2(3,1)))) / coeff2(1,1))), hold on
%      [4:0.02:9], (coeff3(1,1)./(1+10.^(coeff3(2,1)*([4:0.02:9] - coeff3(3,1)))) / coeff3(1,1)), ...
%      [4:0.02:9], (coeff4(1,1)./(1+10.^(coeff4(2,1)*([4:0.02:9] - coeff4(3,1)))) / coeff4(1,1))), hold on

legend('SXB2/3 = 0.5/0.1','SXB2/3 = 1/1',  'location', 'Northwest')%,'SXB2/3 = 0.1/0.5','SXB2/3 = 0.5/0.1',  'location', 'Northwest')
errorbar(pCa_Range, (squeeze(Force(1,1,:)))/ Force(1,1,1), squeeze(Force_std(1,1,:))./Force(1,1,1), 'bo')
errorbar(pCa_Range, (squeeze(Force(2,2,:))/ Force(2,2,1)), squeeze(Force_std(2,2,:))./Force(2,2,1), 'ko')
% errorbar(pCa_Range, squeeze(Force(1,3,:)), squeeze(Force_std(1,3,:)), 'ko')
% errorbar(pCa_Range, squeeze(Force(3,1,:)), squeeze(Force_std(3,1,:)), 'ko')
title('Normalized Force/pCa')
ylabel('Normalized Force (pN)')
xlabel('pCa')
set(gca, 'xdir', 'reverse')
set(gca, 'box', 'off')

