clear all
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No need to mod.
% No Plots, instead creates a table in DataFiles\Kinetic_Rates
% for combinatorial changes in eta, f1, and y0
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pCa = 4.0;

SXB2_Range = [1];
SXB3_Range = [1];
eta_range = [0.6850, 0.822, 0.959];
f1_range = [0.2800, 0.336, 0.392];
y0_range = [20, 35, 60]; 
SlopeL_range = [-100];
SlopeR_range = [20];
atan_max_range = [100];
GausPeak_range = [2000];


prcnt_Full = zeros(length(eta_range),length(f1_range),length(y0_range));  %prcnt_Full @ (SXB2 index, SXB3 index, pCa index)
prcnt_Full_std = zeros(length(eta_range),length(f1_range),length(y0_range));
ton_Full = zeros(length(eta_range),length(f1_range),length(y0_range));
ton_Full_std = zeros(length(eta_range),length(f1_range),length(y0_range));
ton_Back = zeros(length(eta_range),length(f1_range),length(y0_range));
ton_Back_std = zeros(length(eta_range),length(f1_range),length(y0_range));
ton_ms = zeros(length(eta_range),length(f1_range),length(y0_range));
ton_ms_std = zeros(length(eta_range),length(f1_range),length(y0_range));
Force = zeros(length(eta_range),length(f1_range),length(y0_range));
Force_std = zeros(length(eta_range),length(f1_range),length(y0_range));


iSXB2 = 1;
iSXB3 = 1;
ieta = 1;
if1 = 1;
iy0 = 1;
iSlopeL = 1;
iSlopeR = 1;
iatan_max = 1;
iGausPeak = 1;


SXB2 = SXB2_Range(iSXB2);
SXB3 = SXB3_Range(iSXB3);
SlopeL = SlopeL_range(iSlopeL);
SlopeR = SlopeR_range(iSlopeR);
atan_max = atan_max_range(iatan_max);
GausPeak = GausPeak_range(iGausPeak);

filepath = ['J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex\DataFiles\Kinetic_Rates\'];
OutFile=sprintf('%sKinetic_Rates2.txt', filepath);
fid=fopen(OutFile, 'wt' );	%open outfile--tab delimited text
for ieta = 1:length(eta_range)
    eta = eta_range(ieta);
    for if1 = 1:length(f1_range)
        f1 = f1_range(if1);
        for iy0 = 1:length(y0_range)
            y0 = y0_range(iy0);
            
            foldername = ['SXB2=', num2str(SXB2,'%1.2f'),'_SXB3=', num2str(SXB3,'%1.2f'), '_eta=',num2str(eta,'%1.3f'),'_f1=', num2str(f1,'%1.3f'), '_y0=',num2str(y0,'%d'), '_SL=',num2str(SlopeL,'%d'), '_SR=',num2str(SlopeR,'%d'), '_atan=',num2str(atan_max,'%d'), '_GP=',num2str(GausPeak,'%d') filesep];
            pCa = num2str(pCa,'%1.2f');
            filepath = ['J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex\DataFiles\Kinetic_Rates\',foldername,'SS_ton_AVG_pCa_',pCa,'_ROI_A.txt'];
            Data = importdata(filepath);
            prcnt_Full(ieta,if1,iy0) = Data.data(13);
            prcnt_Full_std(ieta,if1,iy0) = Data.data(14);
            ton_Full(ieta,if1,iy0) = Data.data(11);
            ton_Full_std(ieta,if1,iy0) = Data.data(12);
            ton_Back(ieta,if1,iy0) = Data.data(15);
            ton_Back_std(ieta,if1,iy0) = Data.data(16);
            ton_ms(ieta,if1,iy0) = Data.data(7);
            ton_ms_std(ieta,if1,iy0) = Data.data(8);
            
            filepath = ['J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex\DataFiles\Kinetic_Rates\',foldername,'SSData_pCa_',pCa,'_ROI_A.txt'];
            Data = importdata(filepath);
            Force(ieta,if1,iy0) = mean(Data.data(:,3));
            Force_std(ieta,if1,iy0) = std(Data.data(:,3));
            
        end
        DataOut = [y0_range;squeeze(prcnt_Full(ieta,if1,:))';squeeze(prcnt_Full_std(ieta,if1,:))';squeeze(ton_Full(ieta,if1,:))';squeeze(ton_Full_std(ieta,if1,:))';squeeze(ton_Back(ieta,if1,:))';squeeze(ton_Back_std(ieta,if1,:))';squeeze(ton_ms(ieta,if1,:))';squeeze(ton_ms_std(ieta,if1,:))';squeeze(Force(ieta,if1,:))';squeeze(Force_std(ieta,if1,:))']';
        eta_str = num2str(eta_range(ieta));
        f1_str = num2str(f1_range(if1));
        fprintf(fid, ['eta = ' eta_str ' and f1 = ' f1_str '\n']);
        fprintf(fid, 'y0\tprcnt_Full\tPF_std\tton_Full\tTF_std\tton_Back\tTB_std\tton_ms\tt(ms)_std\tForce\tF_std\n');	%header for outfile
        FormatString=[];
        [~, ColOut]=size(DataOut);
        for i=1:ColOut-1 %for all but last
            FormatString=[FormatString, '%4.3f\t'];
        end
        FormatString=[FormatString, '%4.3f\n'];
        fprintf(fid, FormatString, DataOut');
        fprintf(fid, '\n');
        
    end
end
fclose(fid);     %close the file




% % %% Output Excel Summary %%%%%%%%%%%%%%%%%%
% % 
% %eta
% filepath = ['J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex\DataFiles\Kinetic_Rates\'];
% DataOut = [eta_range;squeeze(prcnt_Full(1,1,:,1,1,1,1,1,1))';squeeze(ton_Full(1,1,:,1,1,1,1,1,1))';squeeze(ton_Back(1,1,:,1,1,1,1,1,1))';squeeze(ton_ms(1,1,:,1,1,1,1,1,1))';squeeze(Force(1,1,:,1,1,1,1,1,1))']';
% OutFile=sprintf('%sKinetic_Rates2.txt', filepath);
% % Open File
% fid=fopen(OutFile, 'wt' );	%open outfile--tab delimited text
% %Create the output file column header:
% fprintf(fid, 'eta\tprcnt_Full\tton_Full\tton_Back\tton_ms\tForce\n');	%header for outfile
% FormatString=[];
% [~, ColOut]=size(DataOut);
% for i=1:ColOut-1 %for all but last
%     FormatString=[FormatString, '%4.3f\t'];
% end
% FormatString=[FormatString, '%4.3f\n'];
% fprintf(fid, FormatString, DataOut');
% % %f1
% fprintf(fid, 'f1\tprcnt_Full\tton_Full\tton_Back\tton_ms\tForce\n');	%header for outfile
% DataOut = [f1_range;squeeze(prcnt_Full(1,1,1,:,1,1,1,1,1))';squeeze(ton_Full(1,1,1,:,1,1,1,1,1))';squeeze(ton_Back(1,1,1,:,1,1,1,1,1))';squeeze(ton_ms(1,1,1,:,1,1,1,1,1))';squeeze(Force(1,1,1,:,1,1,1,1,1))']';
% fprintf(fid, FormatString, DataOut');
% %y0
% fprintf(fid, 'y0\tprcnt_Full\tton_Full\tton_Back\tton_ms\tForce\n');	%header for outfile
% DataOut = [y0_range;squeeze(prcnt_Full(1,1,1,1,:,1,1,1,1))';squeeze(ton_Full(1,1,1,1,:,1,1,1,1))';squeeze(ton_Back(1,1,1,1,:,1,1,1,1))';squeeze(ton_ms(1,1,1,1,:,1,1,1,1))';squeeze(Force(1,1,1,1,:,1,1,1,1))']';
% fprintf(fid, FormatString, DataOut');
% %SlopeL
% fprintf(fid, 'SlopeL\tprcnt_Full\tton_Full\tton_Back\tton_ms\tForce\n');	%header for outfile
% DataOut = [SlopeL_range;squeeze(prcnt_Full(1,1,1,1,1,:,1,1,1))';squeeze(ton_Full(1,1,1,1,1,:,1,1,1))';squeeze(ton_Back(1,1,1,1,1,:,1,1,1))';squeeze(ton_ms(1,1,1,1,1,:,1,1,1))';squeeze(Force(1,1,1,1,1,:,1,1,1))']';
% fprintf(fid, FormatString, DataOut');
% %SlopeR
% fprintf(fid, 'SlopeR\tprcnt_Full\tton_Full\tton_Back\tton_ms\tForce\n');	%header for outfile
% DataOut = [SlopeR_range;squeeze(prcnt_Full(1,1,1,1,1,1,:,1,1))';squeeze(ton_Full(1,1,1,1,1,1,:,1,1))';squeeze(ton_Back(1,1,1,1,1,1,:,1,1))';squeeze(ton_ms(1,1,1,1,1,1,:,1,1))';squeeze(Force(1,1,1,1,1,1,:,1,1))']';
% fprintf(fid, FormatString, DataOut');
% %atan_max
% fprintf(fid, 'atan_max\tprcnt_Full\tton_Full\tton_Back\tton_ms\tForce\n');	%header for outfile
% DataOut = [atan_max_range;squeeze(prcnt_Full(1,1,1,1,1,1,1,:,1))';squeeze(ton_Full(1,1,1,1,1,1,1,:,1))';squeeze(ton_Back(1,1,1,1,1,1,1,:,1))';squeeze(ton_ms(1,1,1,1,1,1,1,:,1))';squeeze(Force(1,1,1,1,1,1,1,:,1))']';
% fprintf(fid, FormatString, DataOut');
% %GausPeak
% fprintf(fid, 'GausPeak\tprcnt_Full\tton_Full\tton_Back\tton_ms\tForce\n');	%header for outfile
% DataOut = [GausPeak_range;squeeze(prcnt_Full(1,1,1,1,1,1,1,1,:))';squeeze(ton_Full(1,1,1,1,1,1,1,1,:))';squeeze(ton_Back(1,1,1,1,1,1,1,1,:))';squeeze(ton_ms(1,1,1,1,1,1,1,1,:))';squeeze(Force(1,1,1,1,1,1,1,1,:))']';
% fprintf(fid, FormatString, DataOut');
% fclose(fid);     %close the file


% %% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%
% iSXB2 = 1;
% iSXB3 = 1;
% ieta = 1;
% if1 = 1;
% iy0 = 1;
% iSlopeL = 1;
% iSlopeR = 1;
% iatan_max = 1;
% iGausPeak = 1;
% 
% clf(figure(1))
% subplot(2,2,1)
% scatter(eta_range,prcnt_Full(iSXB2,iSXB3,:,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak))
% xlabel('eta');
% ylabel('prcnt Full');
% ylim([0.1 0.25])
% subplot(2,2,2)
% scatter(eta_range,ton_Full(iSXB2,iSXB3,:,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak))
% xlabel('eta');
% ylabel('ton Full');
% ylim([8 11])
% subplot(2,2,3)
% scatter(eta_range,ton_Back(iSXB2,iSXB3,:,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak))
% xlabel('eta');
% ylabel('ton Back');
% ylim([2 5])
% subplot(2,2,4)
% scatter(eta_range,ton_ms(iSXB2,iSXB3,:,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak))
% xlabel('eta');
% ylabel('ton ms');
% ylim([3 6])
% 
% clf(figure(2))
% subplot(2,2,1)
% scatter(f1_range,prcnt_Full(iSXB2,iSXB3,ieta,:,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak))
% xlabel('f1');
% ylabel('prcnt Full');
% ylim([0.1 0.25])
% subplot(2,2,2)
% scatter(f1_range,ton_Full(iSXB2,iSXB3,ieta,:,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak))
% xlabel('f1');
% ylabel('ton Full');
% ylim([8 11])
% subplot(2,2,3)
% scatter(f1_range,ton_Back(iSXB2,iSXB3,ieta,:,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak))
% xlabel('f1');
% ylabel('ton Back');
% ylim([2 5])
% subplot(2,2,4)
% scatter(f1_range,ton_ms(iSXB2,iSXB3,ieta,:,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak))
% xlabel('f1');
% ylabel('ton ms');
% ylim([3 6])
% 
% clf(figure(3))
% subplot(2,2,1)
% scatter(y0_range,prcnt_Full(iSXB2,iSXB3,ieta,if1,:,iSlopeL,iSlopeR,iatan_max,iGausPeak))
% xlabel('y0');
% ylabel('prcnt Full');
% ylim([0.1 0.25])
% subplot(2,2,2)
% scatter(y0_range,ton_Full(iSXB2,iSXB3,ieta,if1,:,iSlopeL,iSlopeR,iatan_max,iGausPeak))
% xlabel('y0');
% ylabel('ton Full');
% ylim([8 11])
% subplot(2,2,3)
% scatter(y0_range,ton_Back(iSXB2,iSXB3,ieta,if1,:,iSlopeL,iSlopeR,iatan_max,iGausPeak))
% xlabel('y0');
% ylabel('ton Back');
% ylim([2 5])
% subplot(2,2,4)
% scatter(y0_range,ton_ms(iSXB2,iSXB3,ieta,if1,:,iSlopeL,iSlopeR,iatan_max,iGausPeak))
% xlabel('y0');
% ylabel('ton ms');
% ylim([3 6])
% 
% clf(figure(4))
% subplot(2,2,1)
% scatter(SlopeL_range,prcnt_Full(iSXB2,iSXB3,ieta,if1,iy0,:,iSlopeR,iatan_max,iGausPeak))
% xlabel('SlopeL');
% ylabel('prcnt Full');
% ylim([0.1 0.25])
% subplot(2,2,2)
% scatter(SlopeL_range,ton_Full(iSXB2,iSXB3,ieta,if1,iy0,:,iSlopeR,iatan_max,iGausPeak))
% xlabel('SlopeL');
% ylabel('ton Full');
% ylim([8 11])
% subplot(2,2,3)
% scatter(SlopeL_range,ton_Back(iSXB2,iSXB3,ieta,if1,iy0,:,iSlopeR,iatan_max,iGausPeak))
% xlabel('SlopeL');
% ylabel('ton Back');
% ylim([2 5])
% subplot(2,2,4)
% scatter(SlopeL_range,ton_ms(iSXB2,iSXB3,ieta,if1,iy0,:,iSlopeR,iatan_max,iGausPeak))
% xlabel('SlopeL');
% ylabel('ton ms');
% ylim([3 6])
% 
% clf(figure(5))
% subplot(2,2,1)
% scatter(SlopeR_range,prcnt_Full(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,:,iatan_max,iGausPeak))
% xlabel('SlopeR');
% ylabel('prcnt Full');
% ylim([0.1 0.25])
% subplot(2,2,2)
% scatter(SlopeR_range,ton_Full(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,:,iatan_max,iGausPeak))
% xlabel('SlopeR');
% ylabel('ton Full');
% ylim([8 11])
% subplot(2,2,3)
% scatter(SlopeR_range,ton_Back(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,:,iatan_max,iGausPeak))
% xlabel('SlopeR');
% ylabel('ton Back');
% ylim([2 5])
% subplot(2,2,4)
% scatter(SlopeR_range,ton_ms(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,:,iatan_max,iGausPeak))
% xlabel('SlopeR');
% ylabel('ton ms');
% ylim([3 6])
% 
% clf(figure(6))
% subplot(2,2,1)
% scatter(atan_max_range,prcnt_Full(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,:,iGausPeak))
% xlabel('atan max');
% ylabel('prcnt Full');
% ylim([0.1 0.25])
% subplot(2,2,2)
% scatter(atan_max_range,ton_Full(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,:,iGausPeak))
% xlabel('atan max');
% ylabel('ton Full');
% ylim([8 11])
% subplot(2,2,3)
% scatter(atan_max_range,ton_Back(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,:,iGausPeak))
% xlabel('atan max');
% ylabel('ton Back');
% ylim([2 5])
% subplot(2,2,4)
% scatter(atan_max_range,ton_ms(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,:,iGausPeak))
% xlabel('atan max');
% ylabel('ton ms');
% ylim([3 6])
% 
% clf(figure(7))
% subplot(2,2,1)
% scatter(GausPeak_range,prcnt_Full(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,:))
% xlabel('GausPeak');
% ylabel('prcnt Full');
% ylim([0.1 0.25])
% subplot(2,2,2)
% scatter(GausPeak_range,ton_Full(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,:))
% xlabel('GausPeak');
% ylabel('ton Full');
% ylim([8 11])
% subplot(2,2,3)
% scatter(GausPeak_range,ton_Back(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,:))
% xlabel('GausPeak');
% ylabel('ton Back');
% ylim([2 5])
% subplot(2,2,4)
% scatter(GausPeak_range,ton_ms(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,:))
% xlabel('GausPeak');
% ylabel('ton ms');
% ylim([3 6])