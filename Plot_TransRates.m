clear all

%% %%%%%%%%%%% MOD Here %%%%%%%%%%%%%%
new_eta = 0.8;
new_f1 = 0.5;
new_y0 = 35; 
new_SlopeL = -60;
new_SlopeR = 20;
new_atan_max = 50;
new_GausPeak = 2000;
new_S1 = 1;
new_S2 = 1;
new_S3 = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% new_eta = 0.6850;
% new_f1 = 0.2800;
% new_y0 = 20; 
% new_SlopeL = -100; 
% new_SlopeR = 20;
% new_atan_max = 100;
% new_GausPeak = 2000;
% new_S1 = 1;
% new_S2 = 1;
% new_S3 = 1;


%% ORIGINAL VARIABLE PLOT %%%%%
S1 = 1;
S2 = 1;
S3 = 1;

y0=20; % Offset (20)
SlopeL=-100; %Slope of line left of (100)
SlopeR=20; % Slope of line right of 0 (20)
xL=-1; % xposition left of 0;
xR=1; % xposition right of 0;
yL=SlopeL*xL + y0; % yposition left of 0;
yR=SlopeR*xR + y0; % yposition right of 0;
m=(-y0 + yR - (0.5*sqrt( ((xL*(y0-yR)-xR*(y0-yL))/xL)^2 )))/xR;
A=((xL*(y0-yR)-xR*(y0-yL))^2)/(4*(xL^2)*(xR^2));
atan_max = 100;
GausPeak = 2000;

AA=1;   %% k12 offset
BB=A;   %% k20 Left Sloope (3600)
CC=m;   %% k20 Right Slope (-40)
DD=y0;  %% k20 offset? (20)

ATP=5e-3; 		%  Intracellular ATP in Molar
ADP=30e-6;		%  Intracellular ADP in Molar
phos=3e-3;		%  Intracellular Pi in Molar
kxscaler = 3;
J2pNnm=1e21;
kCal2Joule=4.1868e3;
N_Avo=6.022e23;
Gnot=13*kCal2Joule; % Gnot of ATP hydrolysis in Joule/mol
T=273.15+15; % Temperature in K
J2RT=1/(8.314*T); %To convert J/mol to units of RT [R=8.314 J/(mol K) ]
GnotRT=Gnot*J2RT; %Gnot in units of RT ~ 21
dGtot=-GnotRT-log(ATP/(ADP*phos));
RTnm2 = 8.314*T*J2pNnm/N_Avo; %  to convert RT/nm^2 into pN/nm, such that 1 RT/nm^2 =~ 4 pN/nm @ 310K
kRT = kxscaler/RTnm2; % cross brige spring constant in RT/nm^2
colors = ['k'; 'b'; 'r'; 'g'; 'm'; 'y'; 'c'; ];

eta = 0.6850;
f1 = 0.2800;

clf(figure(1))
% for if1 = 1:length(f1_range)
% for ieta = 1:length(eta_range)  
%     eta = eta_range(ieta);
%     f1 = f1_range(if1);
    
    BridgeNRG=abs(eta*dGtot);
    reach = sqrt(BridgeNRG./kRT); % The cross-brige reach, in nm
    
    
    xin=-5:0.01:15;
    for i=1:length(xin)
        G0_N(i) = 0.0;
        G1_N(i) = f1*dGtot + kRT.*((xin(i) - reach).^2);
        G2_N(i) = eta*dGtot + kRT.*(xin(i)^2);
        G0P_N(i) = G0_N(i) + dGtot;
        
        k01(i) = S1*(GausPeak*sqrt(kRT/(2.0*pi)).*exp(-kRT.*((xin(i) - reach).^2)/2.0));
        k10(i) = k01(i).*exp(G1_N(i) - G0_N(i));
        
        k12(i) = S2*(AA + (atan_max./sqrt(kRT)).*(1-tanh(sqrt(kRT).*(xin(i)-reach))));
        
        k21(i) = k12(i).*exp(G2_N(i) - G1_N(i));
        
        k20(i) = S3*(sqrt(kRT).*(sqrt(A*(xin(i)*xin(i))) + m*xin(i)) + y0);
        
        k02(i) = k20(i).*exp((G0_N(i)- G2_N(i)) + dGtot);
    end
    
    subplot(2,2,1)
    hold on;
    plot(xin,G0_N,'k-')
    plot(xin,G1_N,'k-')
    plot(xin,G2_N,'k-')
    plot(xin,G0P_N,'k-')
    axis([-5 15 -50 5])
    xlabel('x(nm)')
    ylabel('Free Energy (RT)')
        
    subplot(2,2,2)
    hold on;
    plot(xin,k01,'k-')
    plot(xin,k10,'k--')
    axis([-5 15 0 1000])
    xlabel('x(nm)')
    ylabel('Transition Rate (s^-^1)')
    
    subplot(2,2,3)
    hold on;
    plot(xin,k12,'k-')
    plot(xin,k21,'k--')
    axis([-5 15 0 1000])
    xlabel('x(nm)')
    ylabel('Transition Rate (s^-^1)')
    
    subplot(2,2,4)
    hold on;
    plot(xin,k20,'k-')
    plot(xin,k02,'k--')
    axis([-5 15 0 1000])
    xlabel('x(nm)')
    ylabel('Transition Rate (s^-^1)')
        
% end
% end



%% New VARIABLE PLOT %%%%%
xL=-1; % xposition left of 0;
xR=1; % xposition right of 0;
yL=new_SlopeL*xL + new_y0; % yposition left of 0;
yR=new_SlopeR*xR + new_y0; % yposition right of 0;
m=(-new_y0 + yR - (0.5*sqrt( ((xL*(new_y0-yR)-xR*(new_y0-yL))/xL)^2 )))/xR;
A=((xL*(new_y0-yR)-xR*(new_y0-yL))^2)/(4*(xL^2)*(xR^2));

AA=1;   %% k12 offset
BB=A;   %% k20 Left Sloope (3600)
CC=m;   %% k20 Right Slope (-40)
DD=new_y0;  %% k20 offset? (20)

ATP=5e-3; 		%  Intracellular ATP in Molar
ADP=30e-6;		%  Intracellular ADP in Molar
phos=3e-3;		%  Intracellular Pi in Molar
kxscaler = 3;
J2pNnm=1e21;
kCal2Joule=4.1868e3;
N_Avo=6.022e23;
Gnot=13*kCal2Joule; % Gnot of ATP hydrolysis in Joule/mol
T=273.15+15; % Temperature in K
J2RT=1/(8.314*T); %To convert J/mol to units of RT [R=8.314 J/(mol K) ]
GnotRT=Gnot*J2RT; %Gnot in units of RT ~ 21
dGtot=-GnotRT-log(ATP/(ADP*phos));
RTnm2 = 8.314*T*J2pNnm/N_Avo; %  to convert RT/nm^2 into pN/nm, such that 1 RT/nm^2 =~ 4 pN/nm @ 310K
kRT = kxscaler/RTnm2; % cross brige spring constant in RT/nm^2
colors = ['k'; 'b'; 'r'; 'g'; 'm'; 'y'; 'c'; ];

% for if1 = 1:length(f1_range)
% for ieta = 1:length(eta_range)  
%     new_eta = eta_range(ieta);
%     new_f1 = f1_range(if1);
    
    BridgeNRG=abs(new_eta*dGtot);
    reach = sqrt(BridgeNRG./kRT); % The cross-brige reach, in nm
    
    
    xin=-5:0.01:15;
    for i=1:length(xin)
        G0_N(i) = 0.0;
        G1_N(i) = new_f1*dGtot + kRT.*((xin(i) - reach).^2);
        G2_N(i) = new_eta*dGtot + kRT.*(xin(i)^2);
        G0P_N(i) = G0_N(i) + dGtot;
        
        k01(i) = new_S1*(new_GausPeak*sqrt(kRT/(2.0*pi)).*exp(-kRT.*((xin(i) - reach).^2)/2.0));
        k10(i) = k01(i).*exp(G1_N(i) - G0_N(i));
        
        k12(i) = new_S2*(AA + (new_atan_max./sqrt(kRT)).*(1-tanh(sqrt(kRT).*(xin(i)-reach))));
        
        k21(i) = k12(i).*exp(G2_N(i) - G1_N(i));
        
        k20(i) = new_S3*(sqrt(kRT).*(sqrt(A*(xin(i)*xin(i))) + m*xin(i)) + new_y0);
        
        k02(i) = k20(i).*exp((G0_N(i)- G2_N(i)) + dGtot);
    end
    
    subplot(2,2,1)
    hold on;
    plot(xin,G0_N,'b-')
    plot(xin,G1_N,'b-')
    plot(xin,G2_N,'b-')
    plot(xin,G0P_N,'b-')
    axis([-5 15 -50 5])
    xlabel('x(nm)')
    ylabel('Free Energy (RT)')
        
    subplot(2,2,2)
    hold on;
    plot(xin,k01,'b-')
    plot(xin,k10,'b--')
    axis([-5 15 0 1000])
    xlabel('x(nm)')
    ylabel('Transition Rate (s^-^1)')
    
    subplot(2,2,3)
    hold on;
    plot(xin,k12,'b-')
    plot(xin,k21,'b--')
    axis([-5 15 0 1000])
    xlabel('x(nm)')
    ylabel('Transition Rate (s^-^1)')
    
    subplot(2,2,4)
    hold on;
    plot(xin,k20,'b-')
    plot(xin,k02,'b--')
    axis([-5 15 0 1000])
    xlabel('x(nm)')
    ylabel('Transition Rate (s^-^1)')
        
% end
% end


xpos = -2;
ypos = 875;
text(xpos,ypos+50,'LEGEND');
text(xpos,ypos,'Original')
text(xpos,ypos-50,'New Value')
temp_legend1 = text(xpos-2,ypos+40,'__');
set(temp_legend1, 'Color', 'k')
set(temp_legend1, 'FontSize', 20)
temp_legend1 = text(xpos-2,ypos-10,'__');
set(temp_legend1, 'Color', 'b')
set(temp_legend1, 'FontSize', 20)




% xpos = -2;
% ypos = 875;
% text(xpos,ypos+50,'LEGEND')
% text(xpos,ypos,'            f1        eta')
% icolor = 1;
% for if1 = 1:length(f1_range)
% for ieta = 1:length(eta_range) 
%     ypos = ypos - 50;
%     text(xpos+2,ypos,[num2str(f1_range(if1),'%1.3f'), '   ', num2str(eta_range(ieta),'%1.3f')]);
%     temp_legend = text(xpos,ypos+40,'__');
%     set(temp_legend, 'Color', colors(icolor))
%     set(temp_legend, 'FontSize', 20)
%     icolor = icolor + 1;
% end
% end