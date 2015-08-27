function []=ton_pCa_Data(InDir, pCa, OutDirectory, SetNumber)

%% Now Gather Our Run Data and Do some Processing
InRun_Data=[InDir 'XBBindingData_pCa_' num2str(pCa, '%3.2f') '.txt'];
disp(['Importing Run Data: ' InRun_Data])
RunData=importdata(InRun_Data);
Runs=RunData.data;

%% Now sort out all of our variables, for how many runs to look at for each
%% set.  Here we have our Run Index (col 1) startin over each time
%% we lanuch a new force trace:
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
        RunsThisTrace=Runs([NewTraceIndex(i,1):1:length(Runs)],:); %% Gather All Runs for this trace
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
        RunsThisTrace=Runs([NewTraceIndex(i,1):1:NewTraceIndex(i+1, 1)-1],:); %% Gather All Runs for this trace
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

OutRunData=[ton_raw, N_raw, ton_full, Pct_full, ton_back, Pct_back, XBF_raw, XBF_full, XBF_back, XBNRG_raw, XBNRG_full, XBNRG_back, Total_XBNRG_raw, Total_XBNRG_full, Total_XBNRG_back];

%% Write out the raw--not averaged--data
% Add path to access the write out the file portion...
OutFile=[OutDirectory 'Set' num2str(SetNumber) '_SS_ton_Raw_pCa_' num2str(pCa, '%3.2f') '.txt']; %% Out Raw Data
disp(['Writing Raw Run Data: ' OutFile])

%% Build Header
head={%'ActF (pN)';...
    %'ATP myo s';...
    %'F per ATP';...
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
fid=fopen(OutFile,'w');	%open outfile--tab delimited text
%Create the output file column header:
fprintf(fid, formatHead);	%header for outfile
% Write Data File and close file
fprintf(fid, FormatString, OutRunData');
fclose(fid);     %close the file