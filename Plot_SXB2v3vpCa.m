clc
clear all
%% MOD HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FOR PLOTTING HILL CURVES
%Trace 1 (Subplot 2 will be normalized to TRACE 2)
T1_SXB2 = 0.5; %0.1, 0.25, 0.5, 1.0, or 1.5
T1_SXB3 = 0.1; %0.1, 0.25, 0.5, 1.0, or 1.5
%Trace 2 
T2_SXB2 = 1.0; %0.1, 0.25, 0.5, 1.0, or 1.5
T2_SXB3 = 1.0; %0.1, 0.25, 0.5, 1.0, or 1.5
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

pCa_Range = [4.0, 4.5, 5.0, 5.5, 5.7, 5.8, 5.9, 6.0, 6.1, 6.25, 6.5, 7.0, 8.0];

SXB2_Range = [0.1 0.25 0.5 1 1.5];
SXB3_Range = [0.1 0.25 0.5 1 1.5];

prcnt_Full = zeros(length(SXB2_Range),length(SXB3_Range),length(pCa_Range));  %prcnt_Full @ (SXB2 index, SXB3 index, pCa index)
ton_Full = zeros(length(SXB2_Range),length(SXB3_Range),length(pCa_Range));
ton_Back = zeros(length(SXB2_Range),length(SXB3_Range),length(pCa_Range));
ton_ms = zeros(length(SXB2_Range),length(SXB3_Range),length(pCa_Range));
Force = zeros(length(SXB2_Range),length(SXB3_Range),length(pCa_Range));
Force_std = zeros(length(SXB2_Range),length(SXB3_Range),length(pCa_Range));
TS_Time = {{},{},{}};  % TimeSeries Time @ (SXB2 index, SXB3 index, pCa index)
TS_Force = {{},{},{}}; % TimeSeries Force @ (SXB2 index, SXB3 index, pCa index)
Ca_level = {{},{},{}};


disp('Working.... This operation may take a few seconds.')

for iSXB2 = 1:length(SXB2_Range)
    SXB2 = num2str(SXB2_Range(iSXB2),'%1.2f');
    for iSXB3 = 1:length(SXB3_Range)
        SXB3 = num2str(SXB3_Range(iSXB3),'%1.2f');
        foldername = ['Muscle=Soleus_Rate=0_HalfSl=1210_SXB2=',SXB2,'_SXB3=',SXB3,filesep];
        for ipCa = 1:length(pCa_Range)
            pCa = num2str(pCa_Range(ipCa),'%1.2f');
            filepath = ['J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex\DataFiles\SXB2v3vpCa\',foldername,'SS_ton_AVG_pCa_',pCa,'_ROI_A.txt'];
            Data = importdata(filepath);
            prcnt_Full(iSXB2,iSXB3,ipCa) = Data.data(13);
            ton_Full(iSXB2,iSXB3,ipCa) = Data.data(11);
            ton_Back(iSXB2,iSXB3,ipCa) = Data.data(15);
            ton_ms(iSXB2,iSXB3,ipCa) = Data.data(7);
            
            filepath = ['J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex\DataFiles\SXB2v3vpCa\',foldername,'SSData_pCa_',pCa,'_ROI_A.txt'];
            Data = importdata(filepath);
            Force(iSXB2,iSXB3,ipCa) = mean(Data.data(:,3));
            Force_std(iSXB2,iSXB3,ipCa) = std(Data.data(:,3));
        end
    end
end

iT1_SXB2 = find(SXB2_Range == T1_SXB2);
iT1_SXB3 = find(SXB3_Range == T1_SXB3);
iT2_SXB2 = find(SXB2_Range == T2_SXB2);
iT2_SXB3 = find(SXB3_Range == T2_SXB3);


clf(figure(1))
subplot(2,2,1)
surf(SXB2_Range,SXB3_Range,prcnt_Full(:,:,1))
title('prcnt Full');
xlabel('SXB2');
ylabel('SXB3');
zlabel('prcnt Full')
subplot(2,2,2)
surf(SXB2_Range,SXB3_Range,ton_ms(:,:,1))
title('ton ms');
xlabel('SXB2');
ylabel('SXB3');
zlabel('ton ms')
subplot(2,2,3)
surf(SXB2_Range,SXB3_Range,ton_Full(:,:,1))
title('ton Full');
xlabel('SXB2');
ylabel('SXB3');
zlabel('ton Full')
subplot(2,2,4)
surf(SXB2_Range,SXB3_Range,ton_Back(:,:,1))
title('ton Back');
xlabel('SXB2');
ylabel('SXB3');
zlabel('ton Back')

clf(figure(2))
subplot(2,2,1)
contour(SXB2_Range,SXB3_Range,prcnt_Full(:,:,1))
title('prcnt Full');
xlabel('SXB2');
ylabel('SXB3');
zlabel('prcnt Full')
subplot(2,2,2)
contour(SXB2_Range,SXB3_Range,ton_ms(:,:,1))
title('ton ms');
xlabel('SXB2');
ylabel('SXB3');
zlabel('ton ms')
subplot(2,2,3)
contour(SXB2_Range,SXB3_Range,ton_Full(:,:,1))
title('ton Full');
xlabel('SXB2');
ylabel('SXB3');
zlabel('ton Full')
subplot(2,2,4)
contour(SXB2_Range,SXB3_Range,ton_Back(:,:,1))
title('ton Back');
xlabel('SXB2');
ylabel('SXB3');
zlabel('ton Back')

% clf(figure(3))
% plot(pCa_Range, squeeze(Force(3,1,:)), pCa_Range, squeeze(Force(4,4,:)))
% % plot(pCa_Range, squeeze(Force(1,3,:))), hold on
% % plot(pCa_Range, squeeze(Force(3,1,:))), hold on
% set(gca,'XDir','Reverse')
% xlabel('pCa');
% ylabel('Force (pN)');
% legend('SXB2/3 = 0.5/0.1','SXB2/3 = 1/1',  'location', 'Northwest')

% 
% Force(:,:,5) = 0; 
% Force_std(:,:,5) = 0;
% pCa_Range(5) = 8.0;

[head1, coeff1]=Process_3ParamHill_Anal_v2(pCa_Range', squeeze(Force(iT1_SXB2,iT1_SXB3,:)));
[head2, coeff2]=Process_3ParamHill_Anal_v2(pCa_Range', squeeze(Force(iT2_SXB2,iT2_SXB3,:)));

% head = {};
% coeff = {};
% for i = 1:length(SXB2_Range)
%     for j = 1:length(SXB3_Range)
%         [head_temp, coeff_temp]=Process_3ParamHill_Anal_v2(pCa_Range', squeeze(Force(iT1_SXB2,iT1_SXB3,:)));
%         head(i,j) = {head_temp};
%         coeff(i,j) = {coeff_temp};
%     end
% end

clf(figure(4))
% subplot(2,2,1)
% for i = 1:length(SXB2_Range)
%     for j = 1:length(SXB3_Range)
%         plot([4:0.02:9], (coeff{i,j}(1,1)./(1+10.^(coeff{i,j}(2,1)*([4:0.02:9] - coeff{i,j}(3,1))))))
%         %errorbar(pCa_Range',squeeze(XB(iSXB2,iSXB3,:,n)),squeeze(XB_std(iSXB2,iSXB3,:,n)),['-',Col{n},'.'],'MarkerSize',15);
%         hold on
%     end
% end


clf(figure(4))
subplot(2,2,1)
plot([4:0.02:9], (coeff1(1,1)./(1+10.^(coeff1(2,1)*([4:0.02:9] - coeff1(3,1))))),...
     [4:0.02:9], (coeff2(1,1)./(1+10.^(coeff2(2,1)*([4:0.02:9] - coeff2(3,1)))))), hold on

legend(sprintf('SXB2/3 = %1.2f/%1.2f',T1_SXB2,T1_SXB3),sprintf('SXB2/3 = %1.2f/%1.2f',T2_SXB2,T2_SXB3),'location', 'Northwest')%,'SXB2/3 = 0.1/0.5','SXB2/3 = 0.5/0.1',  'location', 'Northwest')
errorbar(pCa_Range', squeeze(Force(iT1_SXB2,iT1_SXB3,:)), squeeze(Force_std(iT1_SXB2,iT1_SXB3,:)),'bo')
errorbar(pCa_Range', squeeze(Force(iT2_SXB2,iT2_SXB3,:)), squeeze(Force_std(iT2_SXB2,iT2_SXB3,:)),'ko')

title('Force/pCa')
ylabel('Force (pN)')
xlabel('pCa')
set(gca, 'xdir', 'reverse')
set(gca, 'box', 'off')

h2 = subplot(subplot(2,2,2));
set(gca,'visible','off')
text(0.2,0.4, sprintf('For SXB2 of %1.2f and SXB3 of %1.2f\nFmax = %4.2f (STD %4.2f)\nHill Coeff = %4.2f (STD %4.2f)\npCa50 = %4.2f [%4.2f to %4.2f] (STD %4.2f)\n\nFor SXB2 of %1.2f and SXB3 of %1.2f\nFmax = %4.2f (STD %4.2f)\nHill Coeff = %4.2f (STD %4.2f)\npCa50 = %4.2f [%4.2f to %4.2f] (STD %4.2f)\n',T1_SXB2,T1_SXB3,coeff1(1,1), coeff1(1,2), coeff1(2,1), coeff1(2,2), coeff1(3,1), coeff1(3,4),coeff1(3,5),coeff1(3,2),T2_SXB2,T2_SXB3,coeff2(1,1), coeff2(1,2), coeff2(2,1), coeff2(2,2), coeff2(3,1), coeff2(3,4),coeff2(3,5),coeff2(3,2)), 'parent', h2);


subplot(2,2,3)
plot([4:0.02:9], (coeff1(1,1)./(1+10.^(coeff1(2,1)*([4:0.02:9] - coeff1(3,1)))) / coeff2(1,1)), ...
     [4:0.02:9], (coeff2(1,1)./(1+10.^(coeff2(2,1)*([4:0.02:9] - coeff2(3,1)))) / coeff2(1,1))), hold on
 
legend(sprintf('SXB2/3 = %1.2f/%1.2f',T1_SXB2,T1_SXB3),sprintf('SXB2/3 = %1.2f/%1.2f',T2_SXB2,T2_SXB3),'location', 'Northwest')%,'SXB2/3 = 0.1/0.5','SXB2/3 = 0.5/0.1',  'location', 'Northwest')
errorbar(pCa_Range, squeeze(Force(iT1_SXB2,iT1_SXB3,:))./ max(squeeze(Force(iT2_SXB2,iT2_SXB3,:))), squeeze(Force_std(iT1_SXB2,iT1_SXB3,:))./Force(2,2,1), 'bo')
errorbar(pCa_Range, squeeze(Force(iT2_SXB2,iT2_SXB3,:))./ max(squeeze(Force(iT2_SXB2,iT2_SXB3,:))), squeeze(Force_std(iT2_SXB2,iT2_SXB3,:))./Force(2,2,1), 'ko')
title('Relative Force/pCa')
ylabel('Relative Force (pN)')
xlabel('pCa')
set(gca, 'xdir', 'reverse')
set(gca, 'box', 'off')

subplot(2,2,4)
plot([4:0.02:9], (coeff1(1,1)./(1+10.^(coeff1(2,1)*([4:0.02:9] - coeff1(3,1)))) / coeff1(1,1)), ...
     [4:0.02:9], (coeff2(1,1)./(1+10.^(coeff2(2,1)*([4:0.02:9] - coeff2(3,1)))) / coeff2(1,1))), hold on

legend(sprintf('SXB2/3 = %1.2f/%1.2f',T1_SXB2,T1_SXB3),sprintf('SXB2/3 = %1.2f/%1.2f',T2_SXB2,T2_SXB3),'location', 'Northwest')%,'SXB2/3 = 0.1/0.5','SXB2/3 = 0.5/0.1',  'location', 'Northwest')
errorbar(pCa_Range, squeeze(Force(iT1_SXB2,iT1_SXB3,:))./ max(squeeze(Force(iT1_SXB2,iT1_SXB3,:))), squeeze(Force_std(iT1_SXB2,iT1_SXB3,:))./Force(1,1,1), 'bo')
errorbar(pCa_Range, squeeze(Force(iT2_SXB2,iT2_SXB3,:))./ max(squeeze(Force(iT2_SXB2,iT2_SXB3,:))), squeeze(Force_std(iT2_SXB2,iT2_SXB3,:))./Force(2,2,1), 'ko')
title('Normalized Force/pCa')
ylabel('Normalized Force (pN)')
xlabel('pCa')
set(gca, 'xdir', 'reverse')
set(gca, 'box', 'off')


disp('Done!')