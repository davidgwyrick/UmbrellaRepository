% generates figure 7a and b
function dynamic_twitch_species
addpath('./export_fig','..');
addpath('J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\Land_Niederer\TnImodel-release');
set(0,'defaultlinelinewidth',1.5);
set(0,'DefaultAxesFontSize',10);
set(0,'DefaultTextFontSize',10);

solver = @ode23t;

figure(71); clf; hold on;
for split = {'','_nosplit'}
 splitlbl = char(split);
 ratCai = [0.237 0.236 0.234 0.235 0.245 0.293 0.310 0.363 0.459 0.579 0.703 0.821 0.927 1.019 1.098 1.166 1.222 1.269 1.308 1.340 1.366 1.387 1.402 1.413 1.420 1.423 1.423 1.419 1.413 1.403 1.392 1.378 1.363 1.347 1.329 1.311 1.293 1.274 1.255 1.236 1.218 1.199 1.181 1.163 1.145 1.127 1.110 1.094 1.077 1.061 1.045 1.030 1.014 0.999 0.985 0.970 0.956 0.942 0.928 0.915 0.901 0.888 0.875 0.862 0.850 0.837 0.825 0.813 0.801 0.789 0.778 0.766 0.755 0.744 0.733 0.722 0.712 0.701 0.691 0.681 0.671 0.661 0.651 0.642 0.632 0.623 0.614 0.605 0.596 0.588 0.579 0.571 0.562 0.554 0.547 0.539 0.531 0.524 0.516 0.509 0.502 0.495 0.488 0.482 0.475 0.469 0.462 0.456 0.450 0.444 0.438 0.433 0.427 0.421 0.416 0.411 0.406 0.401 0.396 0.391 0.386 0.381 0.377 0.372 0.368 0.363 0.359 0.355 0.351 0.347 0.343 0.339 0.335 0.332 0.328 0.324 0.321 0.318 0.314 0.311 0.308 0.305 0.301 0.298 0.295 0.292 0.289 0.287 0.284 0.281 0.278 0.276 0.273 0.271 0.268 0.266 0.263 0.261 0.259 0.256 0.254 0.252 0.250 0.248 0.245 0.243 0.241 0.237];
 options = dynamic_xb_prep(['rat' splitlbl]);

 f=@(t,y) make_dy(t,y,interp1(0:167,ratCai,mod(t,167)),options);
 [~,y] = solver(@(t,y) f(t,y),0:1:(6*167), options.y0);
 frat = (options.state_xbj-1)' * y(:,1:options.num_tmxb)';

 mouseCai = 1e-3 * [199.0,198.8,199.9,202.3,212.1,227.5,248.6,279.7,312.8,347.9,378.9,406.1,429.5,445.2,458.3,468.7,476.6,483.8,490.3,494.2,496.3,496.5,494.4,492.1,489.4,487.5,485.9,484.8,483.5,481.6,479.3,475.2,470.7,466.1,461.2,456.7,452.6,449.9,447.0,444.1,440.4,436.1,431.2,426.3,421.8,417.7,415.4,413.9,413.0,411.6,409.1,405.7,402.1,398.5,394.8,392.7,391.1,389.8,388.6,386.9,384.7,380.9,377.1,373.4,369.6,366.4,363.7,361.6,359.2,356.5,352.7,348.9,345.1,341.5,338.5,336.2,334.1,332.2,330.7,328.1,325.2,322.1,318.8,316.0,313.8,312.3,311.0,310.0,308.2,305.9,303.3,300.2,297.4,294.7,293.1,291.9,291.1,289.8,288.3,286.6,283.9,281.4,279.2,277.1,275.7,275.0,274.5,273.6,272.3,269.9,267.5,265.2,263.3,261.7,260.6,259.7,258.8,258.0,256.5,254.9,253.0,250.7,248.7,247.0,245.7,244.9,244.6,243.6,242.3,240.7,238.5,236.5,234.8,233.6,232.7,232.3,231.7,230.8,229.6,227.8,225.8,223.8,222.7,221.8,221.0,220.9,220.4,219.5,218.0,216.5,215.1,213.5,212.3,211.3,210.5,209.9,209.6,208.6,207.5,206.3,204.9,203.6,202.7,201.9,201.5,201.4,200.5 199.0];
  
 options = dynamic_xb_prep(['mouse' splitlbl]);
 f=@(t,y) make_dy(t,y,interp1(0:167,mouseCai,mod(t,167)),options);
 [~,y] = solver(@(t,y) f(t,y),0:1:(6*167), options.y0);
 fmouse = (options.state_xbj-1)' * y(:,1:options.num_tmxb)';

 humanCai=make_ca(0.1399,0.3431,48.2,175.9,343.1); %coppini data
 options = dynamic_xb_prep(['human' splitlbl]);
 f=@(t,y) make_dy(t,y,interp1(0:1000,humanCai,mod(t,1000)),options);
 [~,y] = solver(@(t,y) f(t,y),0:1:(3*1000), options.y0);
 fhuman = (options.state_xbj-1)' * y(:,1:options.num_tmxb)';
     
 K = (120 / (0.25*69));
 if isempty(splitlbl)
   linetype = '-';
 else
   linetype = '--';
 end
 plot(0:500,K*fmouse(end-500:end),linetype,'Color',[0.8 0.2 0.2],'LineWidth',2);
 dx=0;plot(0:500,K*frat(end-dx-500:end-dx),linetype,'Color',[0.2 0.2 0.8],'LineWidth',2);
 plot(0:1000,K*fhuman(end-1000:end),linetype,'Color',[0.2 0.8 0.2],'LineWidth',2);
 
 ylabel('Force (kPa)');
 xlabel('Time (ms)')
end

xlim([0 500]);
legend('q=0.5  Mouse','       Rat','       Human  ','q=1    Mouse','       Rat','       Human  ','Location','NorthEast');
myfig2(figure(71), './Fig7a.png',[3.25 2.5]*2,225);

figure(72); clf; hold on; 
  plot(0:1001,[mouseCai(1:end-1) mouseCai(1:end-1) mouseCai(1:end-1) mouseCai(1:end-1) mouseCai(1:end-1) mouseCai(1:end-1)],'-','Color',[0.8 0.2 0.2],'LineWidth',2);
  plot(0:1001,[ratCai(1:end-1) ratCai(1:end-1) ratCai(1:end-1) ratCai(1:end-1) ratCai(1:end-1) ratCai(1:end-1)],'-','Color',[0.2 0.2 0.8],'LineWidth',2);
  plot(0:1000,humanCai,'-','Color',[0.2 0.8 0.2],'LineWidth',2);
  xlim([0 1000]); legend('Mouse', 'Rat', 'Human');
  xlabel('Time (ms)'); ylabel('[Ca^{2+}]_i \mu{}M');
 
myfig2(figure(72), './Fig7b.png',[3.25 2.5]*2,225);

