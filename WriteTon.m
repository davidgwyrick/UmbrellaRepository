%% Import Force, ATPase, F/ATPase, ton Data
%% Make a ton plot
%% output to outfile:


function []=WriteTon(dt, InDir, pCa, OutDirectory, XB_Fraction, Stats)

%function []=ton_pCa_Data_2013_7_31()

% dt = 0.001;
% pCa = 5.0;
% InDir = 'DataFiles\Muscle=Soleus_Rate=0_HalfSl=1210_SXB2=0.1_SXB3=0.1\';
% OutDirectory = 'DataFiles\Muscle=Soleus_Rate=0_HalfSl=1210_SXB2=0.1_SXB3=0.1\';
% XB_Fraction = 1;


SetNumber = 1;
Figure_Number = 2;

ROIsAZ = 'A':'Z';
ROIRange = ROIsAZ(1:size(Stats.MFMean,3));
% ROIRange = ROIsAZ(1:1);

% for Rate_index = 1:length(RateRange)
%     Rate = RateRange{Rate_index};
%
%     for SL_index = 1:length(SLRange)
%         SL = SLRange{SL_index};
%
%         for Shorten_index = 1:length(ShortenRange)
%             Shorten = ShortenRange{Shorten_index};
%
for ROI_index = 1:length(ROIRange)
    ROI = ROIRange(ROI_index);
    
    %                 for pCa_index = 1:length(pCaRange)
    %                     pCa = pCaRange(pCa_index);
    
    %                     InDir = ['DataFiles\Rate_',Rate,'_SL_',SL,'_',Shorten,'nm\'];
    %                     OutDirectory = ['DataFiles\Rate_',Rate,'_SL_',SL,'_',Shorten,'nm\'];
    
    %% Import Steady-State Force, ATPase, Data
    InSS_Data=[InDir 'TwitchData_pCa_' num2str(pCa, '%3.2f') '_ROI_' ROI '.txt'];
    disp(['Importing Steady State Data: ' InSS_Data])
    SSdata=importdata(InSS_Data);
    SSdata=SSdata.data;
    
    %% Import Data for each particular column
    iActive_Force=5; %% Import the Active force==Think filament force calc
    iATP_per_dt=17;  %% Import ATP per time step
    %% Scale ATP per dt==> ATPase per myosin in simulation
    %% calculate time step
    %dt=Data{1,2}(2,1);
    Raw_N_Bridges=4*3*60; %% (NMYO*N_Thick_Start*XB_Nodes)==>for Vert3Simple Only
    ATP_Scaler=(dt*Raw_N_Bridges*XB_Fraction);
    
    %% Force, ATPase, F/ATPase
    SS_Data_Out=[SSdata(:,iActive_Force), SSdata(:,iATP_per_dt)/ATP_Scaler, SSdata(:,iActive_Force)./(SSdata(:,iATP_per_dt)/ATP_Scaler)];
    
    
    %% Now Gather Our Run Data and Do some Processing
    InRun_Data=[InDir 'XBBindingData_pCa_' num2str(pCa, '%3.2f') '_ROI_' ROI '.txt'];
    disp(['Importing Run Data: ' InRun_Data])
    RunData=importdata(InRun_Data);
    if isstruct(RunData) == 1; % Checking if empty. (If the output of XBBindingData is empty, it won't be input here as a struct, just an array with the headers)
        Runs=RunData.data;
    else
        Runs = NaN(1, 18);
        Runs(1,1) = 1;
    end
    %% Now sort out all of our variables, for how many runs to look at for each
    %% set.  Here we have our Run Index (col 1) startin over each time
    %% we launch a new force trace:
    %%
    NewTraceCol=1; %% column 1
    NewTraceIndex=find(Runs(:, NewTraceCol)==1); %% This will give us a handle
    %% to sort on each Run as a new trace startes for keeping track of ton for
    %% each trace (once it settles out in the end)
    ton_raw=zeros(length(NewTraceIndex), 1);
    ton_full=zeros(length(NewTraceIndex), 1); %% For Full Cycle
    ton_back=zeros(length(NewTraceIndex), 1); %% Did not complete full cycle
    N_raw=zeros(length(NewTraceIndex), 1); % Number of Raw
    Pct_full=zeros(length(NewTraceIndex), 1); % Percent of Full
    Pct_back=zeros(length(NewTraceIndex), 1); %% Percent of backward exit
    % and exited out the back door.
    XBF_raw=zeros(length(NewTraceIndex), 1); %% XB Force (mean)
    XBF_full=zeros(length(NewTraceIndex), 1); %% XB Force Full Cycle
    XBF_back=zeros(length(NewTraceIndex), 1); %% XB Force Back
    XBNRG_raw=zeros(length(NewTraceIndex), 1); %% XB Energy (mean)
    XBNRG_full=zeros(length(NewTraceIndex), 1); %% XB Energy Full Cycle
    XBNRG_back=zeros(length(NewTraceIndex), 1); %% XB Energy Back
    Total_XBNRG_raw=zeros(length(NewTraceIndex), 1); %% Total XB Energy (sum)
    Total_XBNRG_full=zeros(length(NewTraceIndex), 1); %% Total XB Energy Full Cycle
    Total_XBNRG_back=zeros(length(NewTraceIndex), 1); %% Total XB Energy Back
    
    %% Some handles
    i_ton=3;
    iComplete=18;%% The column for completed runs is 18
    ifull=8; % If the 8th column is 1, then it went throgh full cycle
    %% otherwise == 0 for going out of the back end
    iXBstrain=11;
    ikXB=16;
    
    
    
    %% Loop through each trace for a given Ca2+ level
    for i=1:length(NewTraceIndex)
        if i==length(NewTraceIndex)
            if NewTraceIndex(end) == size(Runs,1)
                RunsThisTrace = Runs(NewTraceIndex(i,1),:); %if there is only one row of binding data (Rare, but happens at high pCa)
            else
                RunsThisTrace=Runs([NewTraceIndex(i,1):1:size(Runs,1)],:); %% Gather All Runs for this trace
            end
            index_CompRuns=find(RunsThisTrace(:,iComplete)==1); %% Get the index for all completed runs this trace
            CompRunsThisTrace=RunsThisTrace(index_CompRuns,:); %% Gather All Completed Runs for this trace
            FullRunsThisTrace=CompRunsThisTrace(find(CompRunsThisTrace(:,ifull)==1),:); %% Gather All Completed Runs that finish the cycle
            BackRunsThisTrace=CompRunsThisTrace(find(CompRunsThisTrace(:,ifull)==0),:); %% Gather All Completed Runs that exited backward out the cycle
            ton_raw(i,1)=mean(CompRunsThisTrace(:,i_ton)); %% Gather the time on for all completed runs
            ton_full(i,1)=mean(FullRunsThisTrace(:,i_ton)); %% Gather the time on for all Full runs
            ton_back(i,1)=mean(BackRunsThisTrace(:,i_ton)); %% Gather the time on for all Full runs
            N_raw(i,1)=length(CompRunsThisTrace); % Number of Raw
            Pct_full(i,1)=length(FullRunsThisTrace)/length(CompRunsThisTrace); %% Percent
            Pct_back(i,1)=length(BackRunsThisTrace)/length(CompRunsThisTrace); %% Percent
            XBF_full(i,1)=mean(FullRunsThisTrace(:,iXBstrain).*FullRunsThisTrace(:,ikXB));
            XBF_back(i,1)=mean((BackRunsThisTrace(:,iXBstrain)./BackRunsThisTrace(:,i_ton)).*BackRunsThisTrace(:,ikXB)); %% Some error in accounting if not exit out back--not divided sum yet
            XBF_raw(i,1)=mean( [FullRunsThisTrace(:,iXBstrain).*FullRunsThisTrace(:,ikXB); ...
                (BackRunsThisTrace(:,iXBstrain)./BackRunsThisTrace(:,i_ton)).*BackRunsThisTrace(:,ikXB)] ); %% Some error in accounting if not exit out back--not divided sum yet
            XBNRG_full(i,1)=mean((FullRunsThisTrace(:,iXBstrain).^2).*FullRunsThisTrace(:,ikXB));
            XBNRG_back(i,1)=mean(((BackRunsThisTrace(:,iXBstrain)./BackRunsThisTrace(:,i_ton)).^2).*BackRunsThisTrace(:,ikXB));
            XBNRG_raw(i,1)=mean( [(FullRunsThisTrace(:,iXBstrain).^2).*FullRunsThisTrace(:,ikXB); ...
                ((BackRunsThisTrace(:,iXBstrain)./BackRunsThisTrace(:,i_ton)).^2).*BackRunsThisTrace(:,ikXB)] );
            Total_XBNRG_full(i,1)=sum((FullRunsThisTrace(:,iXBstrain).^2).*FullRunsThisTrace(:,ikXB));
            Total_XBNRG_back(i,1)=sum(((BackRunsThisTrace(:,iXBstrain)./BackRunsThisTrace(:,i_ton)).^2).*BackRunsThisTrace(:,ikXB));
            Total_XBNRG_raw(i,1)=sum( [(FullRunsThisTrace(:,iXBstrain).^2).*FullRunsThisTrace(:,ikXB); ...
                ((BackRunsThisTrace(:,iXBstrain)./BackRunsThisTrace(:,i_ton)).^2).*BackRunsThisTrace(:,ikXB)] );
            
        else
            if (NewTraceIndex(i,1) + 1) == NewTraceIndex(i+1,1)
                RunsThisTrace = Runs(NewTraceIndex(i,1),:); %if there is only one row of binding data (Rare, but happens at high pCa)
            else
                RunsThisTrace=Runs([NewTraceIndex(i,1):1:NewTraceIndex(i+1, 1)-1],:); %% Gather All Runs for this trace
            end
            index_CompRuns=find(RunsThisTrace(:,iComplete)==1); %% Get the index for all completed runs this trace
            CompRunsThisTrace=RunsThisTrace(index_CompRuns,:); %% Gather All Completed Runs for this trace
            FullRunsThisTrace=CompRunsThisTrace(find(CompRunsThisTrace(:,ifull)==1),:); %% Gather All Completed Runs that finish the cycle
            BackRunsThisTrace=CompRunsThisTrace(find(CompRunsThisTrace(:,ifull)==0),:); %% Gather All Completed Runs that exited backward out the cycle
            ton_raw(i,1)=mean(CompRunsThisTrace(:,i_ton)); %% Gather the time on for all completed runs
            ton_full(i,1)=mean(FullRunsThisTrace(:,i_ton)); %% Gather the time on for all Full runs
            ton_back(i,1)=mean(BackRunsThisTrace(:,i_ton)); %% Gather the time on for all Full runs
            N_raw(i,1)=length(CompRunsThisTrace); % Number of Raw
            Pct_full(i,1)=length(FullRunsThisTrace)/length(CompRunsThisTrace); %% Percent
            Pct_back(i,1)=length(BackRunsThisTrace)/length(CompRunsThisTrace); %% Percent
            XBF_full(i,1)=mean(FullRunsThisTrace(:,iXBstrain).*FullRunsThisTrace(:,ikXB));
            XBF_back(i,1)=mean((BackRunsThisTrace(:,iXBstrain)./BackRunsThisTrace(:,i_ton)).*BackRunsThisTrace(:,ikXB)); %% Some error in accounting if not exit out back--not divided sum yet
            XBF_raw(i,1)=mean( [FullRunsThisTrace(:,iXBstrain).*FullRunsThisTrace(:,ikXB); ...
                (BackRunsThisTrace(:,iXBstrain)./BackRunsThisTrace(:,i_ton)).*BackRunsThisTrace(:,ikXB)] ); %% Some error in accounting if not exit out back--not divided sum yet
            XBNRG_full(i,1)=mean((FullRunsThisTrace(:,iXBstrain).^2).*FullRunsThisTrace(:,ikXB));
            XBNRG_back(i,1)=mean(((BackRunsThisTrace(:,iXBstrain)./BackRunsThisTrace(:,i_ton)).^2).*BackRunsThisTrace(:,ikXB));
            XBNRG_raw(i,1)=mean( [(FullRunsThisTrace(:,iXBstrain).^2).*FullRunsThisTrace(:,ikXB); ...
                ((BackRunsThisTrace(:,iXBstrain)./BackRunsThisTrace(:,i_ton)).^2).*BackRunsThisTrace(:,ikXB)] );
            Total_XBNRG_full(i,1)=sum((FullRunsThisTrace(:,iXBstrain).^2).*FullRunsThisTrace(:,ikXB));
            Total_XBNRG_back(i,1)=sum(((BackRunsThisTrace(:,iXBstrain)./BackRunsThisTrace(:,i_ton)).^2).*BackRunsThisTrace(:,ikXB));
            Total_XBNRG_raw(i,1)=sum( [(FullRunsThisTrace(:,iXBstrain).^2).*FullRunsThisTrace(:,ikXB); ...
                ((BackRunsThisTrace(:,iXBstrain)./BackRunsThisTrace(:,i_ton)).^2).*BackRunsThisTrace(:,ikXB)] );
        end
    end
    
    clf(figure(Figure_Number))
    subplot(2,2,1)
    errorbar(1, mean(ton_raw), std(ton_raw), 'ko', 'markerfacecolor', 'k');, hold on
    errorbar(2, mean(ton_full), std(ton_full), 'ko', 'markerfacecolor', 'g');
    errorbar(3, mean(ton_back), std(ton_back), 'ko', 'markerfacecolor', 'm');
    set(gca, 'box', 'off')
    ylabel('t_{on} (ms)')
    xlim([0,4])
    ylim([0,25])
    set(gca, 'xtick', [1,2,3])
    set(gca, 'xticklabel', {'All', 'Full', 'Back Exit'})
    title(['Set: ' num2str(SetNumber) ', pCa=' num2str(pCa)])
    
    subplot(2,2,2)
    errorbar(1, mean(N_raw)/mean(N_raw), std(ton_raw)/mean(N_raw), 'ko', 'markerfacecolor', 'k'); hold on
    errorbar(2, mean(Pct_full), std(Pct_full), 'ko', 'markerfacecolor', 'g');, hold on
    errorbar(3, mean(Pct_back), std(Pct_back), 'ko', 'markerfacecolor', 'm');
    set(gca, 'box', 'off')
    ylabel('Num. Events (percent)')
    xlim([0,4])
    ylim([0,1.2])
    set(gca, 'xtick', [1,2,3])
    set(gca, 'xticklabel', {'All', 'Full', 'Back Exit'})
    title(['N Total/trace ~' num2str(mean(N_raw)) '     ' date])
    
    subplot(2,2,3)
    errorbar(1, mean(XBF_raw), std(XBF_raw), 'ko', 'markerfacecolor', 'k');, hold on
    errorbar(2, mean(XBF_full), std(XBF_full), 'ko', 'markerfacecolor', 'g');
    errorbar(3, mean(XBF_back), std(XBF_back), 'ko', 'markerfacecolor', 'm');
    set(gca, 'box', 'off')
    ylabel('XB Force (pN)')
    xlim([0,4])
    ylim([-10,5])
    set(gca, 'xtick', [1,2,3])
    set(gca, 'xticklabel', {'All', 'Full', 'Back Exit'})
    
    subplot(2,2,4)
    errorbar(1, mean(XBNRG_raw), std(XBNRG_raw), 'ko', 'markerfacecolor', 'k');, hold on
    errorbar(2, mean(XBNRG_full), std(XBNRG_full), 'ko', 'markerfacecolor', 'g');
    errorbar(3, mean(XBNRG_back), std(XBNRG_back), 'ko', 'markerfacecolor', 'm');
    set(gca, 'box', 'off')
    ylabel('XB Energy (pN nm)')
    xlim([0,4])
    ylim([0,100])
    set(gca, 'xtick', [1,2,3])
    set(gca, 'xticklabel', {'All', 'Full', 'Back Exit'})
    
    Figure_Number = Figure_Number + 1;
    
    OutRunData=[ton_raw, N_raw, ton_full, Pct_full, ton_back, Pct_back, XBF_raw, XBF_full, XBF_back, XBNRG_raw, XBNRG_full, XBNRG_back, Total_XBNRG_raw, Total_XBNRG_full, Total_XBNRG_back];
    
    nanarray = nan(length(SS_Data_Out),15);
    nanarray(1:size(OutRunData,1),:) = OutRunData;
    OutRunData = nanarray;
    
    
    %% Write out the raw--not averaged--data
    % Add path to access the write out the file portion...
    OutFile=[OutDirectory 'Twitch_ton_Raw_pCa_' num2str(pCa, '%3.2f') '_ROI_' ROI '.txt']; %% Out Raw Data
    disp(['Writing Raw Run Data: ' OutFile])
    
    %% Build Header
    head={'ActF (pN)';...
        'ATP myo s';...
        'F per ATP';...
        'ton (ms)';...
        'N ton';...
        'ton Full';...
        'prcnt Full';...
        'ton Back';...
        'prcnt Back';...
        'XB_F (pN)';...
        'XB_F Full';...
        'XB_F Back';...
        'XB_En (pN nm)';...
        'XB_En Full';...
        'XB_En Back';...
        'Tot_En (pN nm)';...
        'Tot_En Full';...
        'Tot_En Back';...
        };
    
    % write data file for parameters from fit
    % format the header for writing the output file
    formatHead=[];
    FormatString=[];
    for i=1:length(head)-1
        formatHead=[formatHead, head{i,1}, '\t'];
        FormatString=[FormatString, '%6.4f\t'];
    end
    formatHead=[formatHead, head{end,1}, '\n'];
    FormatString=[FormatString, '%6.4f\n'];
    %% Write Total Force
    fid=fopen(OutFile,'w+');	%open outfile--tab delimited text
    %Create the output file column header:
    fprintf(fid, formatHead);	%header for outfile
    % Write Data File and close file
    fprintf(fid, FormatString, [SS_Data_Out, OutRunData]');
    fclose(fid);     %close the file
    
    %% Write out the Averaged Data
    OutRunData_Avg=[mean(ton_raw), std(ton_raw),...
        nanmean(N_raw), std(N_raw), ...
        nanmean(ton_full), std(ton_full),...
        nanmean(Pct_full), std(Pct_full),...
        nanmean(ton_back), std(ton_back),...
        nanmean(Pct_back), std(Pct_back), ...
        nanmean(XBF_raw), std(XBF_raw),...
        nanmean(XBF_full), std(XBF_full), ...
        nanmean(XBF_back), std(XBF_back), ...
        nanmean(XBNRG_raw), std(XBNRG_raw), ...
        nanmean(XBNRG_full), std(XBNRG_full), ...
        nanmean(XBNRG_back), mean(XBNRG_back), ...
        nanmean(Total_XBNRG_raw), std(Total_XBNRG_raw), ...
        nanmean(Total_XBNRG_full), std(Total_XBNRG_full), ...
        nanmean(Total_XBNRG_back), std(Total_XBNRG_back)];
    
    SS_Data_Out_Avg=[mean(SS_Data_Out(:,1)), std(SS_Data_Out(:,1)),...
        nanmean(SS_Data_Out(:,2)), std(SS_Data_Out(:,2)),...
        nanmean(SS_Data_Out(:,3)), std(SS_Data_Out(:,3))];
    
    
    %% Write out the raw--not averaged--data
    % Add path to access the write out the file portion...
    OutFile2=[OutDirectory 'Twitch_ton_AVG_pCa_' num2str(pCa, '%3.2f') '_ROI_' ROI '.txt']; %% Out Avg Data
    disp(['Writing Raw Run Data: ' OutFile2])
    
    %% Build Header
    head2={'ActF (pN)';...
        'STD ActF';...
        'ATP myo s';...
        'STD ATP';...
        'F per ATP';...
        'STD FpATP';...
        'ton (ms)';...
        'STD ton';...
        'N ton';...
        'STD N';...
        'ton Full';...
        'STD Full';...
        'prcnt Full';...
        'STD pctF';...
        'ton Back';...
        'STD Back';...
        'prcnt Back';...
        'STD pctB';...
        'XB_F (pN)';...
        'STD XB_F';...
        'XB_F Full';...
        'STD XB_FF';...
        'XB_F Back';...
        'STD XB_FB';...
        'XB_En (pN nm)';...
        'STD XB_En';...
        'XB_En Full';...
        'STD XB_EnF';...
        'XB_En Back';...
        'STD XB_EnB';...
        'Tot_En (pN nm)';...
        'STD XBTotEn';...
        'Tot_En Full';...
        'STD XB_TotF';...
        'Tot_En Back';...
        'STD XB_EnE';...
        };
    
    % write data file for parameters from fit
    % format the header for writing the output file
    formatHead=[];
    FormatString=[];
    for i=1:length(head2)-1
        formatHead=[formatHead, head2{i,1}, '\t'];
        FormatString=[FormatString, '%6.4f\t'];
    end
    formatHead=[formatHead, head2{end,1}, '\n'];
    FormatString=[FormatString, '%6.4f\n'];
    %% Write Total Force
    fid=fopen(OutFile2,'w');	%open outfile--tab delimited text
    %Create the output file column header:
    fprintf(fid, formatHead);	%header for outfile
    % Write Data File and close file
    fprintf(fid, FormatString, [SS_Data_Out_Avg, OutRunData_Avg]');
    fclose(fid);     %close the file
    %                 end
end

%         end
%     end
% end
