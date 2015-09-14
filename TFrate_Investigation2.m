%% Plot Specific TF Rate Folders
% David Wyrick
% 9/8/15
% Purpose: To see how changing koff, CaOff, and RuOff affect the Force-pCa
% curves 
%%

Koff_range = [150];
RuOff_range=[100];
CaOff = [100];
pCa_Range = [4.00; 4.50; 5.00; 5.50; 5.70; 5.80; 5.90; 6.00; 6.10; 6.25; 6.50; 7.00];
Pulse_Width_Range = [90]; %number of ms of Ca2+ pulse

% TS_Time=cell(1,9);
% TS_Force=cell(1,9);
% TS_CaLevel=cell(1,9);
% TS_XBfrac=cell(1,9); %XB Fraction Bound
% TS_TF1=cell(1,9); %Actins Ca0
% TS_TF2=cell(1,9); %Actins Ca1
% TS_TF3=cell(1,9); %Actins Ca2

Force=zeros(12,1);
Stdv_Force=zeros(12,1);

%% Loop To gather data
index=1;
for ii=1:length(Koff_range)
    Koff=Koff_range(ii);
%     for jj=1:length(RuOff_range)
       RuOff=RuOff_range;
    for kk=1:length(pCa_Range)
        pCa=pCa_Range(kk);

        foldername=['TFrates_koff=',num2str(Koff),' RuOff=',num2str(RuOff),' CaOff=',num2str(CaOff)];
        filepath=['J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\DataFiles\TF_Rates\',foldername,'\SSData_pCa_',num2str(pCa,'% 10.2f'),'_ROI_A.txt'];
        if pCa == 4,disp(foldername);end
        %disp(pCa)
        Data = importdata(filepath);
        
        Force(kk,ii)=mean(Data.data(:,3));
        Stdv_Force(kk,ii)=std(Data.data(:,3));
        
    end
    end


%% Hill Curves
[head1,coeff1]=Process_3ParamHill_Anal_v2(pCa_Range,Force(:,1));
% [head2,coeff2]=Process_3ParamHill_Anal_v2(pCa_Range,Force(:,2));
% [head3,coeff3]=Process_3ParamHill_Anal_v2(pCa_Range,Force(:,3));

%%
clf(figure(1))
errorbar(pCa_Range,Force(:,1),Stdv_Force(:,1),'-or','LineWidth',.5,'MarkerSize',2),hold on
label1=['koff=150 RuOff=100 CaOff=100 pCa50=',num2str(coeff1(3,1)),' nh=',num2str(coeff1(2,1))];
% errorbar(pCa_Range,Force(:,2),Stdv_Force(:,2),'-og','LineWidth',.5,'MarkerSize',2),hold on
% label2=['koff=25 RuOff=100 CaOff=100 pCa50=',num2str(coeff2(3,1)),' nh=',num2str(coeff2(2,1))];
% errorbar(pCa_Range,Force(:,3),Stdv_Force(:,3),'-ob','LineWidth',.5,'MarkerSize',2),hold on
% label3=['koff=100 RuOff=100 CaOff=100 pCa50=',num2str(coeff3(3,1)),' nh=',num2str(coeff3(2,1))];
%Legend(label1,label2,label3)
legend(label1)
xlabel('pCa')
ylabel('Force (pN)')
axis([3 7 0 900],'on')
set(gca, 'xdir', 'reverse')
title('Force-pCa curve for Different TF Rate values')
fprintf('\nFor koff=150 RuOff=100 CaOff=100\nFmax = %4.2f\nHill Coeff = %4.2f\npCa50 = %4.2f [%4.2f to %4.2f]\n',coeff1(1,1), coeff1(2,1), coeff1(3,1), coeff1(3,4),coeff1(3,5))
% fprintf('\nFor koff=25 RuOff=5 CaOff=100\nFmax = %4.2f\nHill Coeff = %4.2f\npCa50 = %4.2f [%4.2f to %4.2f]\n',coeff2(1,1), coeff2(2,1), coeff2(3,1), coeff2(3,4),coeff2(3,5))
% fprintf('\nFor koff=100 RuOff=5 CaOff=100\nFmax = %4.2f\nHill Coeff = %4.2f\npCa50 = %4.2f [%4.2f to %4.2f]\n',coeff3(1,1), coeff3(2,1), coeff3(3,1), coeff3(3,4),coeff3(3,5))





    
    