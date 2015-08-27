%function Plot_Folders

pCa_Range = [4.0, 4.5, 5.0, 5.5, 5.7, 5.8, 5.9, 6.0, 6.1, 6.25, 6.5, 7.0, 8.0];
%RateRange = [0]; %-

SXB2_Range = [0.5, 1];
SXB3_Range = [0.1, 1];


prcnt_Full = zeros(length(SXB2_Range),length(SXB3_Range),length(pCa_Range));
ton_Full = zeros(length(SXB2_Range),length(SXB3_Range),length(pCa_Range));
ton_Back = zeros(length(SXB2_Range),length(SXB3_Range),length(pCa_Range));
ton_ms = zeros(length(SXB2_Range),length(SXB3_Range),length(pCa_Range));
Force = zeros(length(SXB2_Range),length(SXB3_Range),length(pCa_Range));
Force_std = zeros(length(SXB2_Range),length(SXB3_Range),length(pCa_Range));

%addpath('J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex')

for iSXB2 = 1:length(SXB2_Range)
    SXB2 = num2str(SXB2_Range(iSXB2),'%1.1f');
    for iSXB3 = 1:length(SXB3_Range)
        SXB3 = num2str(SXB3_Range(iSXB3),'%1.1f');
        foldername = ['Muscle=Soleus_Rate=0_HalfSl=1210_SXB2=',SXB2,'_SXB3=',SXB3,'/'];
        for ipCa = 1:length(pCa_Range)
            pCa = num2str(pCa_Range(ipCa),'%1.2f');
            filepath = ['J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex\DataFiles/',foldername,'SS_ton_AVG_pCa_',pCa,'_ROI_A.txt'];
            Data = importdata(filepath);
            prcnt_Full(iSXB2,iSXB3,ipCa) = Data.data(13);
            ton_Full(iSXB2,iSXB3,ipCa) = Data.data(11);
            ton_Back(iSXB2,iSXB3,ipCa) = Data.data(15);
            ton_ms(iSXB2,iSXB3,ipCa) = Data.data(7);
            
            filepath = ['J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex\DataFiles/',foldername,'SSData_pCa_',pCa,'_ROI_A.txt'];
            Data = importdata(filepath);
            Force(iSXB2,iSXB3,ipCa) = mean(Data.data(:,3));
            Force_std(iSXB2,iSXB3,ipCa) = std(Data.data(:,3));
        end
    end
end






% clf(figure(1))
% subplot(2,2,1)
% surf(SXB2_Range,SXB3_Range,prcnt_Full(:,:,1))
% title('prcnt Full');
% xlabel('SXB2');
% ylabel('SXB3');
% zlabel('prcnt Full')
% subplot(2,2,2)
% surf(SXB2_Range,SXB3_Range,ton_ms(:,:,1))
% title('ton ms');
% xlabel('SXB2');
% ylabel('SXB3');
% zlabel('ton ms')
% subplot(2,2,3)
% surf(SXB2_Range,SXB3_Range,ton_Full(:,:,1))
% title('ton Full');
% xlabel('SXB2');
% ylabel('SXB3');
% zlabel('ton Full')
% subplot(2,2,4)
% surf(SXB2_Range,SXB3_Range,ton_Back(:,:,1))
% title('ton Back');
% xlabel('SXB2');
% ylabel('SXB3');
% zlabel('ton Back')
% 
% clf(figure(2))
% subplot(2,2,1)
% contour(SXB2_Range,SXB3_Range,prcnt_Full(:,:,1))
% title('prcnt Full');
% xlabel('SXB2');
% ylabel('SXB3');
% zlabel('prcnt Full')
% subplot(2,2,2)
% contour(SXB2_Range,SXB3_Range,ton_ms(:,:,1))
% title('ton ms');
% xlabel('SXB2');
% ylabel('SXB3');
% zlabel('ton ms')
% subplot(2,2,3)
% contour(SXB2_Range,SXB3_Range,ton_Full(:,:,1))
% title('ton Full');
% xlabel('SXB2');
% ylabel('SXB3');
% zlabel('ton Full')
% subplot(2,2,4)
% contour(SXB2_Range,SXB3_Range,ton_Back(:,:,1))
% title('ton Back');
% xlabel('SXB2');
% ylabel('SXB3');
% zlabel('ton Back')

clf(figure(3))
plot(pCa_Range, squeeze(Force(1,1,:)), pCa_Range, squeeze(Force(2,2,:)))
% plot(pCa_Range, squeeze(Force(1,3,:))), hold on
% plot(pCa_Range, squeeze(Force(3,1,:))), hold on
set(gca,'XDir','Reverse')
xlabel('pCa');
ylabel('Force (pN)');
legend('SXB2/3 = 0.5/0.1','SXB2/3 = 1/1',  'location', 'Northwest')

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

