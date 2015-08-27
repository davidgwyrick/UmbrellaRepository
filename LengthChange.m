%% SETUP SL Length Input

function [SL_EndLength, NSTEPS, ROI_index] = LengthChange(StartLength, NSTEPS, Rate, DataParams)
DataSpan = DataParams.DataSpan;
dt=DataParams.dt; % simulation timestep, in seconds




%type = 1;      %performs a constant length change for the maximal working space possible
                % _________________
                %                  \
                %                   \
                %                    \
                %            |  A | |B|  

type = 2;       %performs a constant length change for a set distance, holds, then reverses to original length and holds again
                % _________________            ________
                %                  \          /
                %                   \        /
                %                    \______/
                %            |  A | |B|  |C| |D|    |E|   %ROI_index


if Rate == 0    %Over-rides any choice made if the rate = 0; no length change is tacked on to the end
    type = 0;
end





switch type
    
    case 0  %No length change and no added steps
        
        SL_EndLength=StartLength*ones(1, NSTEPS);
        
        %Set region for analysis and returns that index
        EndStepStart = floor((1-DataSpan)*NSTEPS);
        ROI_index = {[EndStepStart:NSTEPS]}; % the steps we record more closely because they are steady-state
        
    case 1  %performs a constant length change for the maximal working space possible
        
        WorkingSpace=StartLength - 1119; %StartLength - Thin Length
                
        Slope=(StartLength*(Rate/(dt*1000))); % (dt*1000)=1; Fractional Shortening in ML/second; (dt*1000) only works if dt=1 ms.
        LastShortenSteps= floor(abs((WorkingSpace/(Slope*dt))));
        LengthTransient=[1:LastShortenSteps]*(Slope*dt)+StartLength;
        
        %Put it all together and what do we get
        StartLengthIndex = StartLength*ones(1, NSTEPS);
        SL_EndLength= [StartLengthIndex LengthTransient];
        
        NSTEPS = length(SL_EndLength);
        
        %Set region for analysis and returns that index
        ROI_index = {[floor((length(StartLengthIndex))-(length(StartLengthIndex))*DataSpan):floor((length(StartLengthIndex)))],...     %10 percent steady state
                    [floor((length(StartLengthIndex))+(length(LengthTransient))*0.5):floor((length(StartLengthIndex))+(length(LengthTransient)))]};     %50 percent slope 1
        
        
    case 2  %performs a constant length change for a set distance, holds, then reverses to original length
        
        WorkingSpace = 100; %nm
        HoldTime1 = 250; %ms (steps)
        HoldTime2 = 250; %ms (steps)
                
        %Shorten for set length (WorkingSpace)
        Slope=(StartLength*(Rate/(dt*1000))); % (dt*1000)=1; Fractional Shortening in ML/second; (dt*1000) only works if dt=1 ms.
        FirstShortenSteps= floor(abs((WorkingSpace/(Slope*dt))));
        ShortenIndex = [1:FirstShortenSteps]*(Slope*dt);
        
        %Hold for set time
        HoldValue1 = ShortenIndex(end);
        HoldIndex1 = zeros(1,HoldTime1)+HoldValue1;
        
        %Lengthen for set length (WorkingSpace)
        Slope=(StartLength*(-Rate/(dt*1000))); % (dt*1000)=1; Fractional Shortening in ML/second; (dt*1000) only works if dt=1 ms.
        LastShortenSteps= floor(abs((WorkingSpace/(Slope*dt))));
        LengthenIndex =  ([1:LastShortenSteps]*(Slope*dt))+HoldValue1;
        
        %Hold for set time
        HoldValue2 = LengthenIndex(end);
        HoldIndex2 = zeros(1,HoldTime2)+HoldValue2;
        
        %Put it all together and what do we get
        StartLengthIndex = StartLength*ones(1, NSTEPS);
        LengthTransient = [ShortenIndex HoldIndex1 LengthenIndex HoldIndex2]+StartLength;
        SL_EndLength= [StartLengthIndex LengthTransient];
        
        NSTEPS = length(SL_EndLength);
        
        
        %Set region for analysis and returns that index
        ROI_index = {[floor((length(StartLengthIndex))-(length(StartLengthIndex))*DataSpan):floor((length(StartLengthIndex)))],...     %10 percent steady state
                    [floor((length(StartLengthIndex))+(length(ShortenIndex))*0.5):floor((length(StartLengthIndex))+(length(ShortenIndex)))],...    %50 percent slope 1
                    [floor((length(StartLengthIndex))+(length(ShortenIndex))+(length(HoldIndex1))*0.5):floor((length(StartLengthIndex))+(length(ShortenIndex))+(length(HoldIndex1)))],... %50 percent hold 1
                    [floor((length(StartLengthIndex))+(length(ShortenIndex))+(length(HoldIndex1))+(length(LengthenIndex))*0.5):floor((length(StartLengthIndex))+(length(ShortenIndex))+(length(HoldIndex1))+(length(LengthenIndex)))],... %50 percent slope 2
                    [floor((length(StartLengthIndex))+(length(ShortenIndex))+(length(HoldIndex1))+(length(LengthenIndex))+(length(HoldIndex2))*0.5):floor((length(StartLengthIndex))+(length(ShortenIndex))+(length(HoldIndex1))+(length(LengthenIndex))+(length(HoldIndex2)))]}; %50 percent hold 2
             
        
        %           clf(figure(1))
        %           subplot(2,1,1)
        %           plot(LengthTransient), hold on
        %           subplot(2,1,2)
        %           plot(SL_EndLength)
        
end

