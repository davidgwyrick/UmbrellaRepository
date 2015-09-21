function [CaChange_index] = CaChange(Ca_protocol, New_pCa, NSTEPS, pulse_width, NumTwitch, TimeBtwTwitches,CaProfile_Scalar)

%New_pCa = 5.5;
pulsetrue=0;
numPulseindex=0;
switch Ca_protocol
    %%
    case 'None'
        CaChange_index = (10^(-New_pCa))*ones(NSTEPS,1);
     %%   
    case 'Step'
        CaChange_index = (10^(-4))*ones(NSTEPS,1);
        CaChange_index(1000:end) = 10^(-New_pCa);
    %%    
    case 'Burst'
        CaChange_index = (0.25*10^-7)*ones(NSTEPS,1);
        CaChange_index(100:100+pulse_width) = 10^(-New_pCa);
     %%   
    case 'Train'
        CaChange_index = (0.25*10^-7)*ones(NSTEPS,1);
        CaChange_index(100:100+pulse_width) = 10^(-New_pCa);
        CaChange_index(550:550+pulse_width) = 10^(-New_pCa);
        CaChange_index(1000:1000+pulse_width) = 10^(-New_pCa);
     %%   
    case 'Twitch' %Need to add 'zero value' for Ca Transient; i.e. start at Resting Diastolic [Ca]
        CaChange_index = (0.25*10^-7)*ones(NSTEPS,1);
        Ca_0 = 0.14*10^(-6); % 0.1399*10^(-6)molar Resting Diastolic Calcium concentration 
        Ca_max =10^(-New_pCa); %0.3431*10^(-6)molar is pCa of 6.46 Ca Transient Max concentration at Time to Peak (48.2ms)
        tao=DetermineTaoParam(pulse_width,New_pCa);
        for indy=1:NSTEPS 
%             if indy == 1, CaChange_index(indy)=Ca_0; end
%             if indy ~=1, indy=indy-1; end
            CaChange_index(indy) = Ca_0 + (Ca_max - Ca_0)*((indy-1)/tao)*exp(1-(indy/tao));
        end
%% Still working on the multiple twitch code        
    case 'MultipleTwitch'
        CaChange_index = (0.25*10^-7)*ones(NSTEPS,1);
        CaChange_base = (0.25*10^-7)*ones(NSTEPS,1);
        Ca_0 = 10^(-6); % 0.1399*10^(-6)molar Resting Diastolic Calcium concentration 
        Ca_max =10^(-4); %0.3431*10^(-6)molar is pCa of 6.46 Ca Transient Max concentration at Time to Peak (48.2ms)
        tao=.5*pulse_width;
        for ii=1:NSTEPS,CaChange_base(ii) = Ca_0 + (Ca_max - Ca_0)*(ii/tao)*exp(1-(ii/tao));end
            
        for indy=1:NSTEPS
            if mod((indy-tao),TimeBtwTwitches) == 0 && numPulseindex < NumTwitch && indy-tao~=0, pulsetrue = pulsetrue+1; end
            switch pulsetrue
                case 0
                  CaChange_index(indy) = Ca_0 + (Ca_max - Ca_0)*(indy/tao)*exp(1-(indy/tao));
                  checkin=indy; 
                case 1
                  CaChange_index(indy) = Ca_0 + (Ca_max - Ca_0)*((indy-checkin)/(tao))*exp(1-((indy-checkin)/(tao)));
                  numPulseindex=1; 
                  checkin1=indy; 
                case 2
                  CaChange_index(indy) = Ca_0 + (Ca_max - Ca_0)*((indy-checkin1)/(tao))*exp(1-((indy-checkin1)/(tao)));
                  numPulseindex=2; 
                  checkin2=indy; 
                case 3
                  CaChange_index(indy) = Ca_0 + (Ca_max - Ca_0)*((indy-checkin2)/(tao))*exp(1-((indy-checkin2)/(tao)));
                  numPulseindex=3; 
            end
                      
%             if pulsetrue == 0,
%                 CaChange_index(indy) = Ca_0 + (Ca_max - Ca_0)*(indy/tao)*exp(1-(indy/tao));
%                 checkin=indy;
%             else
%                 CaChange_index(indy) = Ca_0 + (Ca_max - Ca_0)*((indy-checkin)/(tao))*exp(1-((indy-checkin)/(tao)));
%                 numPulseindex=numPulseindex+1;
%             end
        end
%%        
    case 'LandTwitchHuman'
        filename=['HumanCai.txt'];
        filepath=['J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\DataFiles\CaTransient\',filename];
        humanCai=importdata(filepath);
        hmin=min(humanCai);
        CaChange_index = ones(NSTEPS,1);
        CaScalar=10^(-New_pCa)/(max(humanCai)*10^(-6)-hmin*10^(-6));
        
        if CaProfile_Scalar == 1,
            for ii=1:length(humanCai),
                if NSTEPS == 1000 && ii <1001,
                    CaChange_index(ii)=CaScalar*(humanCai(ii)*10^(-6)-hmin*10^(-6));
                	index=ii+1;
                end
            end
            CaChange_index(index:NSTEPS)= 10^(-8);
        else
            for ii=1:length(humanCai),
                if NSTEPS == 1000 && ii <1001,
                    CaChange_index(ii)=humanCai(ii)*10^(-6);
                    index=ii+1;
                end
            end
            CaChange_index(index:NSTEPS)= min(humanCai);
        end
  %%      
    case 'LandTwitchRat'
        filename=['ratCai.txt'];
        filepath=['J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\DataFiles\CaTransient\',filename];
        ratCai=importdata(filepath);
        rmin=min(ratCai);
        CaScalar=10^(-New_pCa)/(max(ratCai)*10^(-6)-rmin*10^(-6));
        CaChange_index = ones(NSTEPS,1);
        
        jj=1;
        if CaProfile_Scalar == 1,
            for ii=1:NSTEPS
                if mod(ii,168)== 0, jj=1;end
                CaChange_index(ii)=CaScalar*(ratCai(jj)*10^(-6)-rmin*10^(-6));
                jj=jj+1;
            end
        else
            for ii=1:NSTEPS
                if mod(ii,168)== 0, jj=1;end
                CaChange_index(ii)=ratCai(jj)*10^(-6);
                jj=jj+1;
            end
        end
        %%
    case 'LandTwitchMouse'
        filename=['mouseCai.txt'];
        filepath=['J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\DataFiles\CaTransient\',filename];
        mouseCai=importdata(filepath);
        mmin=min(mouseCai);
        CaScalar=10^(-New_pCa)/(max(mouseCai)*10^(-6)-mmin*10^(-6));
        CaChange_index = ones(NSTEPS,1);
        
        jj=1;
        if CaProfile_Scalar == 1,
            for ii=1:NSTEPS
                if mod(ii,168)== 0, jj=1;end
                CaChange_index(ii)=CaScalar*(mouseCai(jj)*10^(-6)-mmin*10^(-6)); %moves the Ca Transient down to base Ca concentration of pCa=8 and scales it to desired peak pCa
                jj=jj+1;
            end    
        else
            for ii=1:NSTEPS
                if mod(ii,168)== 0, jj=1;end
                CaChange_index(ii)=mouseCai(jj)*10^(-6); 
                jj=jj+1;
            end   
        end
        
 
        
        
        
end
       
