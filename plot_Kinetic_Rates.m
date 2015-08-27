clear all
%% % MOD HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No modding here. This script outputs both graphs and
% an excel table with the various ton, tfull, etc for changes
% in individual kinetic rates (no combinations).
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pCa = 4.0;

SXB2_Range = [1];
SXB3_Range = [1];
eta_range = [0.6850, 0.822, 0.548];
f1_range = [0.2800, 0.336, 0.224];
y0_range = [20, 35, 60, 100];
SlopeL_range = [-100, -80, -120];
SlopeR_range = [20, 24, 16];
atan_max_range = [100, 80, 120];
GausPeak_range = [2000, 2400, 1600];


prcnt_Full = zeros(length(SXB2_Range),length(SXB3_Range),length(eta_range),length(f1_range),length(y0_range),length(SlopeL_range),length(SlopeR_range),length(atan_max_range),length(GausPeak_range));  %prcnt_Full @ (SXB2 index, SXB3 index, pCa index)
ton_Full = zeros(length(SXB2_Range),length(SXB3_Range),length(eta_range),length(f1_range),length(y0_range),length(SlopeL_range),length(SlopeR_range),length(atan_max_range),length(GausPeak_range));
ton_Back = zeros(length(SXB2_Range),length(SXB3_Range),length(eta_range),length(f1_range),length(y0_range),length(SlopeL_range),length(SlopeR_range),length(atan_max_range),length(GausPeak_range));
ton_ms = zeros(length(SXB2_Range),length(SXB3_Range),length(eta_range),length(f1_range),length(y0_range),length(SlopeL_range),length(SlopeR_range),length(atan_max_range),length(GausPeak_range));
Force = zeros(length(SXB2_Range),length(SXB3_Range),length(eta_range),length(f1_range),length(y0_range),length(SlopeL_range),length(SlopeR_range),length(atan_max_range),length(GausPeak_range));
Force_std = zeros(length(SXB2_Range),length(SXB3_Range),length(eta_range),length(f1_range),length(y0_range),length(SlopeL_range),length(SlopeR_range),length(atan_max_range),length(GausPeak_range));

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


for ieta = 1:length(eta_range)
    eta = eta_range(ieta);
    f1 = f1_range(if1);
    y0 = y0_range(iy0);
    SlopeL = SlopeL_range(iSlopeL);
    SlopeR = SlopeR_range(iSlopeR);
    atan_max = atan_max_range(iatan_max);
    GausPeak = GausPeak_range(iGausPeak);
    
    foldername = ['SXB2=', num2str(SXB2,'%1.2f'),'_SXB3=', num2str(SXB3,'%1.2f'), '_eta=',num2str(eta,'%1.3f'),'_f1=', num2str(f1,'%1.3f'), '_y0=',num2str(y0,'%d'), '_SL=',num2str(SlopeL,'%d'), '_SR=',num2str(SlopeR,'%d'), '_atan=',num2str(atan_max,'%d'), '_GP=',num2str(GausPeak,'%d') filesep];
    pCa = num2str(pCa,'%1.2f');
    filepath = ['J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex\DataFiles\Kinetic_Rates\',foldername,'SS_ton_AVG_pCa_',pCa,'_ROI_A.txt'];
%     filepath = ['DataFiles\Kinetic_Rates\',foldername,'SS_ton_AVG_pCa_',pCa,'_ROI_A.txt'];
    Data = importdata(filepath);
    prcnt_Full(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = Data.data(13);
    ton_Full(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = Data.data(11);
    ton_Back(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = Data.data(15);
    ton_ms(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = Data.data(7);
    
    filepath = ['J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex\DataFiles\Kinetic_Rates\',foldername,'SSData_pCa_',pCa,'_ROI_A.txt'];
%     filepath = ['DataFiles\Kinetic_Rates\',foldername,'SSData_pCa_',pCa,'_ROI_A.txt'];
    Data = importdata(filepath);
    Force(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = mean(Data.data(:,3));
    Force_std(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = std(Data.data(:,3));
end
ieta = 1;
for if1 = 1:length(f1_range)
    eta = eta_range(ieta);
    f1 = f1_range(if1);
    y0 = y0_range(iy0);
    SlopeL = SlopeL_range(iSlopeL);
    SlopeR = SlopeR_range(iSlopeR);
    atan_max = atan_max_range(iatan_max);
    GausPeak = GausPeak_range(iGausPeak);
    
    foldername = ['SXB2=', num2str(SXB2,'%1.2f'),'_SXB3=', num2str(SXB3,'%1.2f'), '_eta=',num2str(eta,'%1.3f'),'_f1=', num2str(f1,'%1.3f'), '_y0=',num2str(y0,'%d'), '_SL=',num2str(SlopeL,'%d'), '_SR=',num2str(SlopeR,'%d'), '_atan=',num2str(atan_max,'%d'), '_GP=',num2str(GausPeak,'%d') filesep];
    pCa = num2str(pCa,'%1.2f');
    filepath = ['J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex\DataFiles\Kinetic_Rates\',foldername,'SS_ton_AVG_pCa_',pCa,'_ROI_A.txt'];
%     filepath = ['DataFiles\Kinetic_Rates\',foldername,'SS_ton_AVG_pCa_',pCa,'_ROI_A.txt'];
    Data = importdata(filepath);
    prcnt_Full(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = Data.data(13);
    ton_Full(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = Data.data(11);
    ton_Back(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = Data.data(15);
    ton_ms(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = Data.data(7);
    
    filepath = ['J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex\DataFiles\Kinetic_Rates\',foldername,'SSData_pCa_',pCa,'_ROI_A.txt'];
%     filepath = ['DataFiles\Kinetic_Rates\',foldername,'SSData_pCa_',pCa,'_ROI_A.txt'];
    Data = importdata(filepath);
    Force(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = mean(Data.data(:,3));
    Force_std(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = std(Data.data(:,3));
end
if1 = 1;
for iy0 = 1:length(y0_range)
    eta = eta_range(ieta);
    f1 = f1_range(if1);
    y0 = y0_range(iy0);
    SlopeL = SlopeL_range(iSlopeL);
    SlopeR = SlopeR_range(iSlopeR);
    atan_max = atan_max_range(iatan_max);
    GausPeak = GausPeak_range(iGausPeak);
    
    foldername = ['SXB2=', num2str(SXB2,'%1.2f'),'_SXB3=', num2str(SXB3,'%1.2f'), '_eta=',num2str(eta,'%1.3f'),'_f1=', num2str(f1,'%1.3f'), '_y0=',num2str(y0,'%d'), '_SL=',num2str(SlopeL,'%d'), '_SR=',num2str(SlopeR,'%d'), '_atan=',num2str(atan_max,'%d'), '_GP=',num2str(GausPeak,'%d') filesep];
    pCa = num2str(pCa,'%1.2f');
    filepath = ['J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex\DataFiles\Kinetic_Rates\',foldername,'SS_ton_AVG_pCa_',pCa,'_ROI_A.txt'];
%     filepath = ['DataFiles\Kinetic_Rates\',foldername,'SS_ton_AVG_pCa_',pCa,'_ROI_A.txt'];
    Data = importdata(filepath);
    prcnt_Full(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = Data.data(13);
    ton_Full(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = Data.data(11);
    ton_Back(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = Data.data(15);
    ton_ms(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = Data.data(7);
    
    filepath = ['J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex\DataFiles\Kinetic_Rates\',foldername,'SSData_pCa_',pCa,'_ROI_A.txt'];
%     filepath = ['DataFiles\Kinetic_Rates\',foldername,'SSData_pCa_',pCa,'_ROI_A.txt'];
    Data = importdata(filepath);
    Force(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = mean(Data.data(:,3));
    Force_std(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = std(Data.data(:,3));
end
iy0 = 1;
for iSlopeL = 1:length(SlopeL_range)
    eta = eta_range(ieta);
    f1 = f1_range(if1);
    y0 = y0_range(iy0);
    SlopeL = SlopeL_range(iSlopeL);
    SlopeR = SlopeR_range(iSlopeR);
    atan_max = atan_max_range(iatan_max);
    GausPeak = GausPeak_range(iGausPeak);
    
    foldername = ['SXB2=', num2str(SXB2,'%1.2f'),'_SXB3=', num2str(SXB3,'%1.2f'), '_eta=',num2str(eta,'%1.3f'),'_f1=', num2str(f1,'%1.3f'), '_y0=',num2str(y0,'%d'), '_SL=',num2str(SlopeL,'%d'), '_SR=',num2str(SlopeR,'%d'), '_atan=',num2str(atan_max,'%d'), '_GP=',num2str(GausPeak,'%d') filesep];
    pCa = num2str(pCa,'%1.2f');
    filepath = ['J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex\DataFiles\Kinetic_Rates\',foldername,'SS_ton_AVG_pCa_',pCa,'_ROI_A.txt'];
%     filepath = ['DataFiles\Kinetic_Rates\',foldername,'SS_ton_AVG_pCa_',pCa,'_ROI_A.txt'];
    Data = importdata(filepath);
    prcnt_Full(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = Data.data(13);
    ton_Full(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = Data.data(11);
    ton_Back(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = Data.data(15);
    ton_ms(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = Data.data(7);
    
    filepath = ['J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex\DataFiles\Kinetic_Rates\',foldername,'SSData_pCa_',pCa,'_ROI_A.txt'];
%     filepath = ['DataFiles\Kinetic_Rates\',foldername,'SSData_pCa_',pCa,'_ROI_A.txt'];
    Data = importdata(filepath);
    Force(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = mean(Data.data(:,3));
    Force_std(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = std(Data.data(:,3));
end
iSlopeL = 1;
for iSlopeR = 1:length(SlopeR_range)
    eta = eta_range(ieta);
    f1 = f1_range(if1);
    y0 = y0_range(iy0);
    SlopeL = SlopeL_range(iSlopeL);
    SlopeR = SlopeR_range(iSlopeR);
    atan_max = atan_max_range(iatan_max);
    GausPeak = GausPeak_range(iGausPeak);
    
    foldername = ['SXB2=', num2str(SXB2,'%1.2f'),'_SXB3=', num2str(SXB3,'%1.2f'), '_eta=',num2str(eta,'%1.3f'),'_f1=', num2str(f1,'%1.3f'), '_y0=',num2str(y0,'%d'), '_SL=',num2str(SlopeL,'%d'), '_SR=',num2str(SlopeR,'%d'), '_atan=',num2str(atan_max,'%d'), '_GP=',num2str(GausPeak,'%d') filesep];
    pCa = num2str(pCa,'%1.2f');
    filepath = ['J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex\DataFiles\Kinetic_Rates\',foldername,'SS_ton_AVG_pCa_',pCa,'_ROI_A.txt'];
%     filepath = ['DataFiles\Kinetic_Rates\',foldername,'SS_ton_AVG_pCa_',pCa,'_ROI_A.txt'];
    Data = importdata(filepath);
    prcnt_Full(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = Data.data(13);
    ton_Full(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = Data.data(11);
    ton_Back(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = Data.data(15);
    ton_ms(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = Data.data(7);
    
    filepath = ['J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex\DataFiles\Kinetic_Rates\',foldername,'SSData_pCa_',pCa,'_ROI_A.txt'];
%     filepath = ['DataFiles\Kinetic_Rates\',foldername,'SSData_pCa_',pCa,'_ROI_A.txt'];
    Data = importdata(filepath);
    Force(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = mean(Data.data(:,3));
    Force_std(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = std(Data.data(:,3));
end
iSlopeR = 1;
for iatan_max = 1:length(atan_max_range)
    eta = eta_range(ieta);
    f1 = f1_range(if1);
    y0 = y0_range(iy0);
    SlopeL = SlopeL_range(iSlopeL);
    SlopeR = SlopeR_range(iSlopeR);
    atan_max = atan_max_range(iatan_max);
    GausPeak = GausPeak_range(iGausPeak);
    
    foldername = ['SXB2=', num2str(SXB2,'%1.2f'),'_SXB3=', num2str(SXB3,'%1.2f'), '_eta=',num2str(eta,'%1.3f'),'_f1=', num2str(f1,'%1.3f'), '_y0=',num2str(y0,'%d'), '_SL=',num2str(SlopeL,'%d'), '_SR=',num2str(SlopeR,'%d'), '_atan=',num2str(atan_max,'%d'), '_GP=',num2str(GausPeak,'%d') filesep];
    pCa = num2str(pCa,'%1.2f');
    filepath = ['J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex\DataFiles\Kinetic_Rates\',foldername,'SS_ton_AVG_pCa_',pCa,'_ROI_A.txt'];
%     filepath = ['DataFiles\Kinetic_Rates\',foldername,'SS_ton_AVG_pCa_',pCa,'_ROI_A.txt'];
    Data = importdata(filepath);
    prcnt_Full(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = Data.data(13);
    ton_Full(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = Data.data(11);
    ton_Back(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = Data.data(15);
    ton_ms(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = Data.data(7);
    
    filepath = ['J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex\DataFiles\Kinetic_Rates\',foldername,'SSData_pCa_',pCa,'_ROI_A.txt'];
%     filepath = ['DataFiles\Kinetic_Rates\',foldername,'SSData_pCa_',pCa,'_ROI_A.txt'];
    Data = importdata(filepath);
    Force(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = mean(Data.data(:,3));
    Force_std(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = std(Data.data(:,3));
end
iatan_max = 1;
for iGausPeak = 1:length(GausPeak_range)
    eta = eta_range(ieta);
    f1 = f1_range(if1);
    y0 = y0_range(iy0);
    SlopeL = SlopeL_range(iSlopeL);
    SlopeR = SlopeR_range(iSlopeR);
    atan_max = atan_max_range(iatan_max);
    GausPeak = GausPeak_range(iGausPeak);
    
    foldername = ['SXB2=', num2str(SXB2,'%1.2f'),'_SXB3=', num2str(SXB3,'%1.2f'), '_eta=',num2str(eta,'%1.3f'),'_f1=', num2str(f1,'%1.3f'), '_y0=',num2str(y0,'%d'), '_SL=',num2str(SlopeL,'%d'), '_SR=',num2str(SlopeR,'%d'), '_atan=',num2str(atan_max,'%d'), '_GP=',num2str(GausPeak,'%d') filesep];
    pCa = num2str(pCa,'%1.2f');
    filepath = ['J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex\DataFiles\Kinetic_Rates\',foldername,'SS_ton_AVG_pCa_',pCa,'_ROI_A.txt'];
%     filepath = ['DataFiles\Kinetic_Rates\',foldername,'SS_ton_AVG_pCa_',pCa,'_ROI_A.txt'];
    Data = importdata(filepath);
    prcnt_Full(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = Data.data(13);
    ton_Full(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = Data.data(11);
    ton_Back(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = Data.data(15);
    ton_ms(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = Data.data(7);
    
    filepath = ['J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex\DataFiles\Kinetic_Rates\',foldername,'SSData_pCa_',pCa,'_ROI_A.txt'];
%     filepath = ['DataFiles\Kinetic_Rates\',foldername,'SSData_pCa_',pCa,'_ROI_A.txt'];
    Data = importdata(filepath);
    Force(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = mean(Data.data(:,3));
    Force_std(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak) = std(Data.data(:,3));
end
iGausPeak = 1;

%% Output Excel Summary %%%%%%%%%%%%%%%%%%

%eta
filepath = ['J:\TannerBertGroup\GraduateStudents\Axel_Fenwick\Sims\Vert_3Alex\DataFiles\Kinetic_Rates\'];
DataOut = [eta_range;squeeze(prcnt_Full(1,1,:,1,1,1,1,1,1))';squeeze(ton_Full(1,1,:,1,1,1,1,1,1))';squeeze(ton_Back(1,1,:,1,1,1,1,1,1))';squeeze(ton_ms(1,1,:,1,1,1,1,1,1))';squeeze(Force(1,1,:,1,1,1,1,1,1))']';
OutFile=sprintf('%sKinetic_Rates.txt', filepath);
% Open File
fid=fopen(OutFile, 'wt' );	%open outfile--tab delimited text
%Create the output file column header:
fprintf(fid, 'eta\tprcnt_Full\tton_Full\tton_Back\tton_ms\tForce\n');	%header for outfile
FormatString=[];
[~, ColOut]=size(DataOut);
for i=1:ColOut-1 %for all but last
    FormatString=[FormatString, '%4.3f\t'];
end
FormatString=[FormatString, '%4.3f\n'];
fprintf(fid, FormatString, DataOut');
%f1
fprintf(fid, 'f1\tprcnt_Full\tton_Full\tton_Back\tton_ms\tForce\n');	%header for outfile
DataOut = [f1_range;squeeze(prcnt_Full(1,1,1,:,1,1,1,1,1))';squeeze(ton_Full(1,1,1,:,1,1,1,1,1))';squeeze(ton_Back(1,1,1,:,1,1,1,1,1))';squeeze(ton_ms(1,1,1,:,1,1,1,1,1))';squeeze(Force(1,1,1,:,1,1,1,1,1))']';
fprintf(fid, FormatString, DataOut');
%y0
fprintf(fid, 'y0\tprcnt_Full\tton_Full\tton_Back\tton_ms\tForce\n');	%header for outfile
DataOut = [y0_range;squeeze(prcnt_Full(1,1,1,1,:,1,1,1,1))';squeeze(ton_Full(1,1,1,1,:,1,1,1,1))';squeeze(ton_Back(1,1,1,1,:,1,1,1,1))';squeeze(ton_ms(1,1,1,1,:,1,1,1,1))';squeeze(Force(1,1,1,1,:,1,1,1,1))']';
fprintf(fid, FormatString, DataOut');
%SlopeL
fprintf(fid, 'SlopeL\tprcnt_Full\tton_Full\tton_Back\tton_ms\tForce\n');	%header for outfile
DataOut = [SlopeL_range;squeeze(prcnt_Full(1,1,1,1,1,:,1,1,1))';squeeze(ton_Full(1,1,1,1,1,:,1,1,1))';squeeze(ton_Back(1,1,1,1,1,:,1,1,1))';squeeze(ton_ms(1,1,1,1,1,:,1,1,1))';squeeze(Force(1,1,1,1,1,:,1,1,1))']';
fprintf(fid, FormatString, DataOut');
%SlopeR
fprintf(fid, 'SlopeR\tprcnt_Full\tton_Full\tton_Back\tton_ms\tForce\n');	%header for outfile
DataOut = [SlopeR_range;squeeze(prcnt_Full(1,1,1,1,1,1,:,1,1))';squeeze(ton_Full(1,1,1,1,1,1,:,1,1))';squeeze(ton_Back(1,1,1,1,1,1,:,1,1))';squeeze(ton_ms(1,1,1,1,1,1,:,1,1))';squeeze(Force(1,1,1,1,1,1,:,1,1))']';
fprintf(fid, FormatString, DataOut');
%atan_max
fprintf(fid, 'atan_max\tprcnt_Full\tton_Full\tton_Back\tton_ms\tForce\n');	%header for outfile
DataOut = [atan_max_range;squeeze(prcnt_Full(1,1,1,1,1,1,1,:,1))';squeeze(ton_Full(1,1,1,1,1,1,1,:,1))';squeeze(ton_Back(1,1,1,1,1,1,1,:,1))';squeeze(ton_ms(1,1,1,1,1,1,1,:,1))';squeeze(Force(1,1,1,1,1,1,1,:,1))']';
fprintf(fid, FormatString, DataOut');
%GausPeak
fprintf(fid, 'GausPeak\tprcnt_Full\tton_Full\tton_Back\tton_ms\tForce\n');	%header for outfile
DataOut = [GausPeak_range;squeeze(prcnt_Full(1,1,1,1,1,1,1,1,:))';squeeze(ton_Full(1,1,1,1,1,1,1,1,:))';squeeze(ton_Back(1,1,1,1,1,1,1,1,:))';squeeze(ton_ms(1,1,1,1,1,1,1,1,:))';squeeze(Force(1,1,1,1,1,1,1,1,:))']';
fprintf(fid, FormatString, DataOut');
fclose(fid);     %close the file


%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%
iSXB2 = 1;
iSXB3 = 1;
ieta = 1;
if1 = 1;
iy0 = 1;
iSlopeL = 1;
iSlopeR = 1;
iatan_max = 1;
iGausPeak = 1;

clf(figure(1))
subplot(2,2,1)
scatter(eta_range,prcnt_Full(iSXB2,iSXB3,:,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak))
xlabel('eta');
ylabel('prcnt Full');
ylim([0.1 0.25])
subplot(2,2,2)
scatter(eta_range,ton_Full(iSXB2,iSXB3,:,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak))
xlabel('eta');
ylabel('ton Full');
ylim([8 11])
subplot(2,2,3)
scatter(eta_range,ton_Back(iSXB2,iSXB3,:,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak))
xlabel('eta');
ylabel('ton Back');
ylim([2 5])
subplot(2,2,4)
scatter(eta_range,ton_ms(iSXB2,iSXB3,:,if1,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak))
xlabel('eta');
ylabel('ton ms');
ylim([3 6])

clf(figure(2))
subplot(2,2,1)
scatter(f1_range,prcnt_Full(iSXB2,iSXB3,ieta,:,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak))
xlabel('f1');
ylabel('prcnt Full');
ylim([0.1 0.25])
subplot(2,2,2)
scatter(f1_range,ton_Full(iSXB2,iSXB3,ieta,:,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak))
xlabel('f1');
ylabel('ton Full');
ylim([8 11])
subplot(2,2,3)
scatter(f1_range,ton_Back(iSXB2,iSXB3,ieta,:,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak))
xlabel('f1');
ylabel('ton Back');
ylim([2 5])
subplot(2,2,4)
scatter(f1_range,ton_ms(iSXB2,iSXB3,ieta,:,iy0,iSlopeL,iSlopeR,iatan_max,iGausPeak))
xlabel('f1');
ylabel('ton ms');
ylim([3 6])

clf(figure(3))
subplot(2,2,1)
scatter(y0_range,prcnt_Full(iSXB2,iSXB3,ieta,if1,:,iSlopeL,iSlopeR,iatan_max,iGausPeak))
xlabel('y0');
ylabel('prcnt Full');
ylim([0.1 0.25])
subplot(2,2,2)
scatter(y0_range,ton_Full(iSXB2,iSXB3,ieta,if1,:,iSlopeL,iSlopeR,iatan_max,iGausPeak))
xlabel('y0');
ylabel('ton Full');
ylim([8 11])
subplot(2,2,3)
scatter(y0_range,ton_Back(iSXB2,iSXB3,ieta,if1,:,iSlopeL,iSlopeR,iatan_max,iGausPeak))
xlabel('y0');
ylabel('ton Back');
ylim([2 5])
subplot(2,2,4)
scatter(y0_range,ton_ms(iSXB2,iSXB3,ieta,if1,:,iSlopeL,iSlopeR,iatan_max,iGausPeak))
xlabel('y0');
ylabel('ton ms');
ylim([3 6])

clf(figure(4))
subplot(2,2,1)
scatter(SlopeL_range,prcnt_Full(iSXB2,iSXB3,ieta,if1,iy0,:,iSlopeR,iatan_max,iGausPeak))
xlabel('SlopeL');
ylabel('prcnt Full');
ylim([0.1 0.25])
subplot(2,2,2)
scatter(SlopeL_range,ton_Full(iSXB2,iSXB3,ieta,if1,iy0,:,iSlopeR,iatan_max,iGausPeak))
xlabel('SlopeL');
ylabel('ton Full');
ylim([8 11])
subplot(2,2,3)
scatter(SlopeL_range,ton_Back(iSXB2,iSXB3,ieta,if1,iy0,:,iSlopeR,iatan_max,iGausPeak))
xlabel('SlopeL');
ylabel('ton Back');
ylim([2 5])
subplot(2,2,4)
scatter(SlopeL_range,ton_ms(iSXB2,iSXB3,ieta,if1,iy0,:,iSlopeR,iatan_max,iGausPeak))
xlabel('SlopeL');
ylabel('ton ms');
ylim([3 6])

clf(figure(5))
subplot(2,2,1)
scatter(SlopeR_range,prcnt_Full(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,:,iatan_max,iGausPeak))
xlabel('SlopeR');
ylabel('prcnt Full');
ylim([0.1 0.25])
subplot(2,2,2)
scatter(SlopeR_range,ton_Full(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,:,iatan_max,iGausPeak))
xlabel('SlopeR');
ylabel('ton Full');
ylim([8 11])
subplot(2,2,3)
scatter(SlopeR_range,ton_Back(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,:,iatan_max,iGausPeak))
xlabel('SlopeR');
ylabel('ton Back');
ylim([2 5])
subplot(2,2,4)
scatter(SlopeR_range,ton_ms(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,:,iatan_max,iGausPeak))
xlabel('SlopeR');
ylabel('ton ms');
ylim([3 6])

clf(figure(6))
subplot(2,2,1)
scatter(atan_max_range,prcnt_Full(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,:,iGausPeak))
xlabel('atan max');
ylabel('prcnt Full');
ylim([0.1 0.25])
subplot(2,2,2)
scatter(atan_max_range,ton_Full(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,:,iGausPeak))
xlabel('atan max');
ylabel('ton Full');
ylim([8 11])
subplot(2,2,3)
scatter(atan_max_range,ton_Back(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,:,iGausPeak))
xlabel('atan max');
ylabel('ton Back');
ylim([2 5])
subplot(2,2,4)
scatter(atan_max_range,ton_ms(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,:,iGausPeak))
xlabel('atan max');
ylabel('ton ms');
ylim([3 6])

clf(figure(7))
subplot(2,2,1)
scatter(GausPeak_range,prcnt_Full(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,:))
xlabel('GausPeak');
ylabel('prcnt Full');
ylim([0.1 0.25])
subplot(2,2,2)
scatter(GausPeak_range,ton_Full(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,:))
xlabel('GausPeak');
ylabel('ton Full');
ylim([8 11])
subplot(2,2,3)
scatter(GausPeak_range,ton_Back(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,:))
xlabel('GausPeak');
ylabel('ton Back');
ylim([2 5])
subplot(2,2,4)
scatter(GausPeak_range,ton_ms(iSXB2,iSXB3,ieta,if1,iy0,iSlopeL,iSlopeR,iatan_max,:))
xlabel('GausPeak');
ylabel('ton ms');
ylim([3 6])