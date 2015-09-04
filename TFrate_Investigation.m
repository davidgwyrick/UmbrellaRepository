%% Plot TF Rate Folders
% David Wyrick
% 9/4/15
% Purpose: To see how changing koff, CaOff, and RuOff affect the timeseries
% data of a 90ms twitch at pCa of 4.0
%%

Koff_range = [5 25 100];
RuOff_range = [5 25 100];
CaOff_range = [5 25 100];
pCaV = [4];
Pulse_Width_Range = [90]; %number of ms of Ca2+ pulse

TS_Time=cell(1,27);
TS_Force=cell(1,27);
TS_CaLevel=cell(1,27);
TS_XBfrac=cell(1,27); %XB Fraction Bound
TS_TF1=cell(1,27); %Actins Ca0
TS_TF2=cell(1,27); %Actins Ca1
TS_TF3=cell(1,27); %Actins Ca2

%% Loop To gather data
index=1;
for ii=1:length(Koff_range)
    Koff=Koff_range(ii);
    for jj=1:length(RuOff_range)
        RuOff=RuOff_range(jj);
        for kk=1:length(CaOff_range)
            CaOff=CaOff_range(kk);
            foldername=['TFrates_koff=',num2str(Koff),' RuOff=',num2str(RuOff),' CaOff=',num2str(CaOff)];
            filepath=['J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\DataFiles\TF_Rates\',foldername,'\TimeSeriesAvg_pCa_4.00.txt'];
            disp(foldername)
            Data = importdata(filepath);
            TS_Time(index)={Data.data(:,1)};
            TS_Force(index)={Data.data(:,3)};
            TS_CaLevel(index)={Data.data(:,12)};
            TS_XBfrac(index)={Data.data(:,5)};
            TS_TF1(index)={Data.data(:,8)}; 
            TS_TF2(index)={Data.data(:,9)};
            TS_TF3(index)={Data.data(:,10)};
            index=index+1;
        end
    end
end

%% Find Max values and Time to 50%/90% Relaxation

TF_Relax=zeros(27,2);
MaxForce=zeros(27,2);
for ii=1:27
    Force=TS_Force{ii};
    Time=TS_Time{ii};
    RT90_found=0;
    RT50_found=0;
    RT90_ind=0;
    RT50_ind=0;
    for jj=1:length(Force),
        [mm mm_index] =max(Force);
        mnew=mean(Force(mm_index-3:mm_index+3));
        MaxForce(ii,2)=mnew;
        MaxForce(ii,1)=mm;
        if jj > mm_index && (Force(jj)/mnew)< .50 && RT50_found == 0; 
            RT50_ind=jj;
            RT50_found=1;
            TF_Relax(ii,1)= Time(RT50_ind);
        end
        if jj > mm_index && (Force(jj)/mnew)< .10 && RT90_found == 0; 
            RT90_ind=jj;
            RT90_found=1;
            TF_Relax(ii,2)= Time(RT90_ind);
        end
               
    end
end

%% Plotting koff=5 
clf(figure(2))
subplot(3,1,1);
hold on;
plot(TS_Time{1},TS_Force{1},'-m','LineWidth',.25,'MarkerSize',2),hold on
label1=['koff=5 RuOff=5 CaOff=5'];
plot(TS_Time{2},TS_Force{2},'--c','LineWidth',.25,'MarkerSize',2),hold on
label2=['koff=5 RuOff=5 CaOff=25'];
plot(TS_Time{3},TS_Force{3},':r','LineWidth',.25,'MarkerSize',2),hold on
label3=['koff=5 RuOff=5 CaOff=100'];
Legend(label1,label2,label3)
xlabel('Time(s)')
ylabel('Force (pN)')
xlim([0 1])
title('Constant koff=5')


subplot(3,1,2);
hold on;
plot(TS_Time{4},TS_Force{4},'-.g','LineWidth',.25,'MarkerSize',2),hold on
label4=['koff=5 RuOff=25 CaOff 5'];
plot(TS_Time{5},TS_Force{5},'-b','LineWidth',.25,'MarkerSize',2),hold on
label5=['koff=5 RuOff=25 CaOff=25'];
plot(TS_Time{6},TS_Force{6},'--k','LineWidth',.25,'MarkerSize',2),hold on
label6=['koff=5 RuOff=25 CaOff=100'];
legend(label4,label5,label6)
xlabel('Time(s)')
ylabel('Force (pN)')
xlim([0 1])

subplot(3,1,3);
hold on;
plot(TS_Time{7},TS_Force{7},':k','LineWidth',.25,'MarkerSize',2),hold on
label7=['koff=5 RuOff=100 CaOff=5'];
plot(TS_Time{8},TS_Force{8},'-.r','LineWidth',.25,'MarkerSize',2),hold on
label8=['koff=5 RuOff=100 CaOff=25'];
plot(TS_Time{9},TS_Force{9},'--b','LineWidth',.25,'MarkerSize',2),hold on
label9=['koff=5 RuOff=100 CaOff=100'];
legend(label7,label8,label9)
xlabel('Time(s)')
ylabel('Force (pN)')
xlim([0 1])

%% Plotting koff=25

clf(figure(3))
subplot(3,1,1);
hold on;
plot(TS_Time{10},TS_Force{10},'-m','LineWidth',.25,'MarkerSize',2),hold on
label1=['koff=25 RuOff=5 CaOff=5'];
plot(TS_Time{11},TS_Force{11},'--c','LineWidth',.25,'MarkerSize',2),hold on
label2=['koff=25 RuOff=5 CaOff=25'];
plot(TS_Time{12},TS_Force{12},':r','LineWidth',.25,'MarkerSize',2),hold on
label3=['koff=25 RuOff=5 CaOff=100'];
Legend(label1,label2,label3)
xlabel('Time(s)')
ylabel('Force (pN)')
xlim([0 1])
title('Constant koff=25')


subplot(3,1,2);
hold on;
plot(TS_Time{13},TS_Force{13},'-.g','LineWidth',.25,'MarkerSize',2),hold on
label4=['koff=25 RuOff=25 CaOff 5'];
plot(TS_Time{14},TS_Force{14},'-b','LineWidth',.25,'MarkerSize',2),hold on
label5=['koff=25 RuOff=25 CaOff=25'];
plot(TS_Time{15},TS_Force{15},'--k','LineWidth',.25,'MarkerSize',2),hold on
label6=['koff=25 RuOff=25 CaOff=100'];
legend(label4,label5,label6)
xlabel('Time(s)')
ylabel('Force (pN)')
xlim([0 1])

subplot(3,1,3);
hold on;
plot(TS_Time{16},TS_Force{16},':k','LineWidth',.25,'MarkerSize',2),hold on
label7=['koff=25 RuOff=100 CaOff=5'];
plot(TS_Time{17},TS_Force{17},'-.r','LineWidth',.25,'MarkerSize',2),hold on
label8=['koff=25 RuOff=100 CaOff=25'];
plot(TS_Time{18},TS_Force{18},'--b','LineWidth',.25,'MarkerSize',2),hold on
label9=['koff=25 RuOff=100 CaOff=100'];
legend(label7,label8,label9)
xlabel('Time(s)')
ylabel('Force (pN)')
xlim([0 1])

%% Plotting koff=100

clf(figure(4))
subplot(3,1,1);
hold on;
plot(TS_Time{19},TS_Force{19},'-m','LineWidth',.25,'MarkerSize',2),hold on
label1=['koff=100 RuOff=5 CaOff=5'];
plot(TS_Time{20},TS_Force{20},'--c','LineWidth',.25,'MarkerSize',2),hold on
label2=['koff=100 RuOff=5 CaOff=25'];
plot(TS_Time{21},TS_Force{21},':r','LineWidth',.25,'MarkerSize',2),hold on
label3=['koff=100 RuOff=5 CaOff=100'];
Legend(label1,label2,label3)
xlabel('Time(s)')
ylabel('Force (pN)')
xlim([0 1])
title('Constant koff=100')


subplot(3,1,2);
hold on;
plot(TS_Time{22},TS_Force{22},'-.g','LineWidth',.25,'MarkerSize',2),hold on
label4=['koff=100 RuOff=25 CaOff 5'];
plot(TS_Time{23},TS_Force{23},'-b','LineWidth',.25,'MarkerSize',2),hold on
label5=['koff=100 RuOff=25 CaOff=25'];
plot(TS_Time{24},TS_Force{24},'--k','LineWidth',.25,'MarkerSize',2),hold on
label6=['koff=100 RuOff=25 CaOff=100'];
legend(label4,label5,label6)
xlabel('Time(s)')
ylabel('Force (pN)')
xlim([0 1])

subplot(3,1,3);
hold on;
plot(TS_Time{25},TS_Force{25},':k','LineWidth',.25,'MarkerSize',2),hold on
label7=['koff=100 RuOff=100 CaOff=5'];
plot(TS_Time{26},TS_Force{26},'-.r','LineWidth',.25,'MarkerSize',2),hold on
label8=['koff=100 RuOff=100 CaOff=25'];
plot(TS_Time{27},TS_Force{27},'--b','LineWidth',.25,'MarkerSize',2),hold on
label9=['koff=100 RuOff=100 CaOff=100'];
legend(label7,label8,label9)
xlabel('Time(s)')
ylabel('Force (pN)')
xlim([0 1])
%% RuOff=5
% clf(figure(12))
% subplot(3,1,1);
% hold on;
% plot(TS_Time{1},TS_Force{1},'-m','LineWidth',.25,'MarkerSize',2),hold on
% label1=['koff=5 RuOff=5 CaOff=5'];
% plot(TS_Time{2},TS_Force{2},'--c','LineWidth',.25,'MarkerSize',2),hold on
% label2=['koff=5 RuOff=5 CaOff=25'];
% plot(TS_Time{3},TS_Force{3},':r','LineWidth',.25,'MarkerSize',2),hold on
% label3=['koff=5 RuOff=5 CaOff=100'];
% Legend(label1,label2,label3)
% xlabel('Time(s)')
% ylabel('Force (pN)')
% xlim([0 1])
% title('Constant koff=5 and RuOff=5')
% 
% subplot(3,1,2);
% hold on;
% plot(TS_Time{10},TS_Force{10},'-.g','LineWidth',.25,'MarkerSize',2),hold on
% label11=['koff=25 RuOff=5 CaOff=5'];
% plot(TS_Time{11},TS_Force{11},'-b','LineWidth',.25,'MarkerSize',2),hold on
% label12=['koff=25 RuOff=5 CaOff=25'];
% plot(TS_Time{12},TS_Force{12},'--k','LineWidth',.25,'MarkerSize',2),hold on
% label13=['koff=25 RuOff=5 CaOff=100'];
% Legend(label11,label12,label13)
% xlabel('Time(s)')
% ylabel('Force (pN)')
% xlim([0 1])
% title('Constant koff=25 and RuOff=5')
% 
% subplot(3,1,3);
% hold on;
% plot(TS_Time{19},TS_Force{19},'-m','LineWidth',.25,'MarkerSize',2),hold on
% label21=['koff=100 RuOff=5 CaOff=5'];
% plot(TS_Time{20},TS_Force{20},'--c','LineWidth',.25,'MarkerSize',2),hold on
% label22=['koff=100 RuOff=5 CaOff=25'];
% plot(TS_Time{21},TS_Force{21},':r','LineWidth',.25,'MarkerSize',2),hold on
% label23=['koff=100 RuOff=5 CaOff=100'];
% Legend(label21,label22,label23)
% xlabel('Time(s)')
% ylabel('Force (pN)')
% xlim([0 1])
% title('Constant koff=100 and RuOff=5')
% 
% %% RuOff = 25
% clf(figure(20))
% subplot(3,1,1);
% hold on;
% plot(TS_Time{4},TS_Force{4},'-m','LineWidth',.25,'MarkerSize',2),hold on
% label4=['koff=5 RuOff=25 CaOff 5'];
% plot(TS_Time{5},TS_Force{5},'--c','LineWidth',.25,'MarkerSize',2),hold on
% label5=['koff=5 RuOff=25 CaOff=25'];
% plot(TS_Time{6},TS_Force{6},':r','LineWidth',.25,'MarkerSize',2),hold on
% label6=['koff=5 RuOff=25 CaOff=100'];
% legend(label4,label5,label6)
% xlabel('Time(s)')
% ylabel('Force (pN)')
% xlim([0 1])
% title('Constant RuOff=25 and koff=5')
% subplot(3,1,2);
% hold on;
% plot(TS_Time{13},TS_Force{13},'-.g','LineWidth',.25,'MarkerSize',2),hold on
% label4=['koff=25 RuOff=25 CaOff 5'];
% plot(TS_Time{14},TS_Force{14},'-b','LineWidth',.25,'MarkerSize',2),hold on
% label5=['koff=25 RuOff=25 CaOff=25'];
% plot(TS_Time{15},TS_Force{15},'--k','LineWidth',.25,'MarkerSize',2),hold on
% label6=['koff=25 RuOff=25 CaOff=100'];
% legend(label4,label5,label6)
% xlabel('Time(s)')
% ylabel('Force (pN)')
% xlim([0 1])
% title('Constant RuOff=25 and koff=25')
% subplot(3,1,3);
% hold on;
% plot(TS_Time{22},TS_Force{22},':g','LineWidth',.25,'MarkerSize',2),hold on
% label4=['koff=100 RuOff=25 CaOff 5'];
% plot(TS_Time{23},TS_Force{23},'-.r','LineWidth',.25,'MarkerSize',2),hold on
% label5=['koff=100 RuOff=25 CaOff=25'];
% plot(TS_Time{24},TS_Force{24},'--b','LineWidth',.25,'MarkerSize',2),hold on
% label6=['koff=100 RuOff=25 CaOff=100'];
% legend(label4,label5,label6)
% xlabel('Time(s)')
% ylabel('Force (pN)')
% xlim([0 1])
% title('Constant RuOff=25 and koff=100')
% 
% %% RuOff=100
% clf(figure(21))
% subplot(3,1,1);
% hold on;
% plot(TS_Time{7},TS_Force{7},'-m','LineWidth',.25,'MarkerSize',2),hold on
% label4=['koff=5 RuOff=100 CaOff 5'];
% plot(TS_Time{8},TS_Force{8},'--c','LineWidth',.25,'MarkerSize',2),hold on
% label5=['koff=5 RuOff=100 CaOff=25'];
% plot(TS_Time{9},TS_Force{9},':r','LineWidth',.25,'MarkerSize',2),hold on
% label6=['koff=5 RuOff=100 CaOff=100'];
% legend(label4,label5,label6)
% xlabel('Time(s)')
% ylabel('Force (pN)')
% xlim([0 1])
% title('Constant RuOff=100 and koff=5')
% 
% subplot(3,1,2);
% hold on;
% plot(TS_Time{16},TS_Force{16},'-.g','LineWidth',.25,'MarkerSize',2),hold on
% label4=['koff=25 RuOff=100 CaOff 5'];
% plot(TS_Time{17},TS_Force{17},'-b','LineWidth',.25,'MarkerSize',2),hold on
% label5=['koff=25 RuOff=100 CaOff=25'];
% plot(TS_Time{18},TS_Force{18},'--k','LineWidth',.25,'MarkerSize',2),hold on
% label6=['koff=25 RuOff=100 CaOff=100'];
% legend(label4,label5,label6)
% xlabel('Time(s)')
% ylabel('Force (pN)')
% xlim([0 1])
% title('Constant RuOff=100 and koff=25')
% 
% subplot(3,1,3);
% hold on;
% plot(TS_Time{25},TS_Force{25},':g','LineWidth',.25,'MarkerSize',2),hold on
% label4=['koff=100 RuOff=100 CaOff 5'];
% plot(TS_Time{26},TS_Force{26},'-.r','LineWidth',.25,'MarkerSize',2),hold on
% label5=['koff=100 RuOff=100 CaOff=25'];
% plot(TS_Time{27},TS_Force{27},'--b','LineWidth',.25,'MarkerSize',2),hold on
% label6=['koff=100 RuOff=100 CaOff=100'];
% legend(label4,label5,label6)
% xlabel('Time(s)')
% ylabel('Force (pN)')
% xlim([0 1])
% title('Constant RuOff=100 and koff=100')
% 
% %% CaOff=5
% clf(figure(5))
% subplot(3,1,1);
% hold on;
% plot(TS_Time{1},TS_Force{1},'-m','LineWidth',.25,'MarkerSize',2),hold on
% label1=['koff=5 RuOff=5 CaOff=5'];
% plot(TS_Time{4},TS_Force{4},'--c','LineWidth',.25,'MarkerSize',2),hold on
% label2=['koff=5 RuOff=25 CaOff=5'];
% plot(TS_Time{7},TS_Force{7},':r','LineWidth',.25,'MarkerSize',2),hold on
% label3=['koff=5 RuOff=100 CaOff=5'];
% Legend(label1,label2,label3)
% xlabel('Time(s)')
% ylabel('Force (pN)')
% xlim([0 1])
% title('Constant CaOff=5 and koff=5')
% 
% 
% subplot(3,1,2);
% hold on;
% plot(TS_Time{10},TS_Force{10},'-.g','LineWidth',.25,'MarkerSize',2),hold on
% label4=['koff=25 RuOff=5 CaOff=5'];
% plot(TS_Time{13},TS_Force{13},'-b','LineWidth',.25,'MarkerSize',2),hold on
% label5=['koff=25 RuOff=25 CaOff=5'];
% plot(TS_Time{16},TS_Force{16},'--k','LineWidth',.25,'MarkerSize',2),hold on
% label6=['koff=25 RuOff=100 CaOff=5'];
% legend(label4,label5,label6)
% xlabel('Time(s)')
% ylabel('Force (pN)')
% title('Constant CaOff=5 and koff=25')
% xlim([0 1])
% 
% subplot(3,1,3);
% hold on;
% plot(TS_Time{19},TS_Force{19},':g','LineWidth',.25,'MarkerSize',2),hold on
% label7=['koff=100 RuOff=5 CaOff=5'];
% plot(TS_Time{22},TS_Force{22},'-.r','LineWidth',.25,'MarkerSize',2),hold on
% label8=['koff=100 RuOff=25 CaOff=5'];
% plot(TS_Time{25},TS_Force{25},'--b','LineWidth',.25,'MarkerSize',2),hold on
% label9=['koff=100 RuOff=100 CaOff=5'];
% legend(label7,label8,label9)
% xlabel('Time(s)')
% ylabel('Force (pN)')
% title('Constant CaOff=5 and koff=100')
% xlim([0 1])
% 
% %% Plotting CaOff=25
% 
% clf(figure(6))
% subplot(3,1,1);
% hold on;
% plot(TS_Time{2},TS_Force{2},'-m','LineWidth',.25,'MarkerSize',2),hold on
% label1=['koff=5 RuOff=5 CaOff=25'];
% plot(TS_Time{5},TS_Force{5},'--c','LineWidth',.25,'MarkerSize',2),hold on
% label2=['koff=5 RuOff=25 CaOff=25'];
% plot(TS_Time{8},TS_Force{8},':r','LineWidth',.25,'MarkerSize',2),hold on
% label3=['koff=5 RuOff=100 CaOff=25'];
% Legend(label1,label2,label3)
% xlabel('Time(s)')
% ylabel('Force (pN)')
% xlim([0 1])
% title('Constant CaOff=25 and koff=5')
% 
% 
% subplot(3,1,2);
% hold on;
% plot(TS_Time{11},TS_Force{11},'-.g','LineWidth',.25,'MarkerSize',2),hold on
% label4=['koff=25 RuOff=5 CaOff=25'];
% plot(TS_Time{14},TS_Force{14},'-b','LineWidth',.25,'MarkerSize',2),hold on
% label5=['koff=25 RuOff=25 CaOff=25'];
% plot(TS_Time{17},TS_Force{17},'--k','LineWidth',.25,'MarkerSize',2),hold on
% label6=['koff=25 RuOff=100 CaOff=25'];
% legend(label4,label5,label6)
% xlabel('Time(s)')
% ylabel('Force (pN)')
% title('Constant CaOff=25 and koff=25')
% xlim([0 1])
% 
% subplot(3,1,3);
% hold on;
% plot(TS_Time{20},TS_Force{20},':g','LineWidth',.25,'MarkerSize',2),hold on
% label7=['koff=100 RuOff=5 CaOff=25'];
% plot(TS_Time{23},TS_Force{23},'-.r','LineWidth',.25,'MarkerSize',2),hold on
% label8=['koff=100 RuOff=25 CaOff=25'];
% plot(TS_Time{26},TS_Force{26},'--b','LineWidth',.25,'MarkerSize',2),hold on
% label9=['koff=100 RuOff=100 CaOff=25'];
% legend(label7,label8,label9)
% xlabel('Time(s)')
% ylabel('Force (pN)')
% title('Constant CaOff=25 and koff=100')
% xlim([0 1])
% 
% %% Plotting CaOff=100
% 
% clf(figure(7))
% subplot(3,1,1);
% hold on;
% plot(TS_Time{3},TS_Force{3},'-m','LineWidth',.25,'MarkerSize',2),hold on
% label1=['koff=5 RuOff=5 CaOff=100'];
% plot(TS_Time{6},TS_Force{6},'--c','LineWidth',.25,'MarkerSize',2),hold on
% label2=['koff=5 RuOff=25 CaOff=100'];
% plot(TS_Time{9},TS_Force{9},':r','LineWidth',.25,'MarkerSize',2),hold on
% label3=['koff=5 RuOff=100 CaOff=100'];
% Legend(label1,label2,label3)
% xlabel('Time(s)')
% ylabel('Force (pN)')
% xlim([0 1])
% title('Constant CaOff=100 and koff=5')
% 
% 
% subplot(3,1,2);
% hold on;
% plot(TS_Time{12},TS_Force{12},'-.g','LineWidth',.25,'MarkerSize',2),hold on
% label4=['koff=25 RuOff=5 CaOff=100'];
% plot(TS_Time{15},TS_Force{15},'-b','LineWidth',.25,'MarkerSize',2),hold on
% label5=['koff=25 RuOff=25 CaOff=100'];
% plot(TS_Time{18},TS_Force{18},'--k','LineWidth',.25,'MarkerSize',2),hold on
% label6=['koff=25 RuOff=100 CaOff=100'];
% legend(label4,label5,label6)
% xlabel('Time(s)')
% ylabel('Force (pN)')
% title('Constant CaOff=100 and koff=25')
% xlim([0 1])
% 
% subplot(3,1,3);
% hold on;
% plot(TS_Time{21},TS_Force{21},':g','LineWidth',.25,'MarkerSize',2),hold on
% label7=['koff=100 RuOff=5 CaOff=100'];
% plot(TS_Time{24},TS_Force{24},'-.r','LineWidth',.25,'MarkerSize',2),hold on
% label8=['koff=100 RuOff=25 CaOff=100'];
% plot(TS_Time{27},TS_Force{27},'--b','LineWidth',.25,'MarkerSize',2),hold on
% label9=['koff=100 RuOff=100 CaOff=100'];
% legend(label7,label8,label9)
% xlabel('Time(s)')
% ylabel('Force (pN)')
% title('Constant CaOff=25 and koff=100')
% xlim([0 1])

%% Master Plot
% plot(TS_Time{1},TS_Force{1},'-m','LineWidth',.25,'MarkerSize',2),hold on
% label1=['koff=5 RuOff=5 CaOff=5'];
% plot(TS_Time{2},TS_Force{2},'--c','LineWidth',.25,'MarkerSize',2),hold on
% label2=['koff=5 RuOff=5 CaOff=25'];
% plot(TS_Time{3},TS_Force{3},':r','LineWidth',.25,'MarkerSize',2),hold on
% label3=['koff=5 RuOff=5 CaOff=100'];
% plot(TS_Time{13},TS_Force{13},'-.g','LineWidth',.25,'MarkerSize',2),hold on
% label4=['koff=25 RuOff=25 CaOff 5'];
% plot(TS_Time{14},TS_Force{14},'-b','LineWidth',.25,'MarkerSize',2),hold on
% label5=['koff=25 RuOff=25 CaOff=25'];
% plot(TS_Time{15},TS_Force{15},'--k','LineWidth',.25,'MarkerSize',2),hold on
% label6=['koff=25 RuOff=25 CaOff=100'];
% plot(TS_Time{7},TS_Force{7},':y','LineWidth',.25,'MarkerSize',2),hold on
% label7=['koff=5 RuOff=100 CaOff=5'];
% plot(TS_Time{8},TS_Force{8},'-.r','LineWidth',.25,'MarkerSize',2),hold on
% label8=['koff=5 RuOff=100 CaOff=25'];
% plot(TS_Time{9},TS_Force{9},'--b','LineWidth',.25,'MarkerSize',2),hold on
% label9=['koff=5 RuOff=100 CaOff=100'];
% 
% plot(TS_Time{10},TS_Force{10},'+m','LineWidth',.25,'MarkerSize',2),hold on
% label11=['koff=25 RuOff=5 CaOff=5'];
% plot(TS_Time{11},TS_Force{11},'oc','LineWidth',.25,'MarkerSize',2),hold on
% label12=['koff=25 RuOff=5 CaOff=25'];
% plot(TS_Time{12},TS_Force{12},'*r','LineWidth',.25,'MarkerSize',2),hold on
% label13=['koff=25 RuOff=5 CaOff=100'];
% plot(TS_Time{13},TS_Force{13},'.g','LineWidth',.25,'MarkerSize',2),hold on
% label14=['koff=25 RuOff=25 CaOff 5'];
% plot(TS_Time{14},TS_Force{14},'xb','LineWidth',.25,'MarkerSize',2),hold on
% label15=['koff=25 RuOff=25 CaOff=25'];
% plot(TS_Time{15},TS_Force{15},'sk','LineWidth',.25,'MarkerSize',2),hold on
% label16=['koff=25 RuOff=25 CaOff=100'];
% plot(TS_Time{16},TS_Force{16},'dy','LineWidth',.25,'MarkerSize',2),hold on
% label17=['koff=25 RuOff=100 CaOff=5'];
% plot(TS_Time{17},TS_Force{17},'^r','LineWidth',.25,'MarkerSize',2),hold on
% label18=['koff=25 RuOff=100 CaOff=25'];
% plot(TS_Time{18},TS_Force{18},'pb','LineWidth',.25,'MarkerSize',2),hold on
% label19=['koff=25 RuOff=100 CaOff=100'];
% 
% plot(TS_Time{19},TS_Force{19},'-r','LineWidth',.25,'MarkerSize',2),hold on
% label21=['koff=100 RuOff=5 CaOff=5'];
% plot(TS_Time{20},TS_Force{20},'-b','LineWidth',.25,'MarkerSize',2),hold on
% label22=['koff=100 RuOff=5 CaOff=25'];
% plot(TS_Time{21},TS_Force{21},'--y','LineWidth',.25,'MarkerSize',2),hold on
% label23=['koff=100 RuOff=5 CaOff=100'];
% plot(TS_Time{22},TS_Force{22},'--g','LineWidth',.25,'MarkerSize',2),hold on
% label24=['koff=100 RuOff=25 CaOff 5'];
% plot(TS_Time{23},TS_Force{23},'-b','LineWidth',.25,'MarkerSize',2),hold on
% label25=['koff=100 RuOff=25 CaOff=25'];
% plot(TS_Time{24},TS_Force{24},'>g','LineWidth',.25,'MarkerSize',2),hold on
% label26=['koff=100 RuOff=25 CaOff=100'];
% plot(TS_Time{25},TS_Force{25},'<y','LineWidth',.25,'MarkerSize',2),hold on
% label27=['koff=100 RuOff=100 CaOff=5'];
% plot(TS_Time{26},TS_Force{26},'-.r','LineWidth',.25,'MarkerSize',2),hold on
% label28=['koff=100 RuOff=100 CaOff=25'];
% plot(TS_Time{27},TS_Force{27},'hb','LineWidth',.25,'MarkerSize',2),hold on
% label29=['koff=100 RuOff=100 CaOff=100'];
% 
% legend(label1,label2,label3,label4,label5,label6,label7,label8,label9,label11,label12,label13,label14,label15,label16,label17,label18,label19,label21,label22,label23,label24,label25,label26,label27,label28,label29)
% xlabel('Time(s)')
% ylabel('Force (pN)')
% xlim([0 1])