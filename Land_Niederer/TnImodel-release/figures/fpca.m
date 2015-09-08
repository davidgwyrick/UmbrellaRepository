% generates figure 6
function fpca
addpath('./export_fig','..');
set(0,'defaultlinelinewidth',1.5);
set(0,'DefaultAxesFontSize',10);
set(0,'DefaultTextFontSize',10);


options = dynamic_xb_prep;
SE = options.SE;

% fpca with xb inhibitor/enhancer
options.pca   = 2:0.01:8;
options.cas =  10.^-options.pca(end:-1:1) * 1e6;
    
KDM_factor = [1 3 0.75]; 
KDMS =  options.KDM * KDM_factor;

f = zeros(length(options.cas),length(KDMS));
nf = f;
    
nC_red  = repmat((0:26)',1,70);
nxb_red = repmat((0:69),27,1);
   
for i=1:length(KDMS)
  KDM = KDMS(i);
  p_red = @(ca) SE .* (1 ./ options.KDA.^(options.n-nC_red)) .* (1 + ca/options.KDC).^(options.n-nC_red) .* (1 + ca/options.KDC + ca/(options.KDC*options.KDI)).^nC_red .* (1 ./ KDM.^nxb_red);
  for ca_ix = 1:length(options.cas)
    p = p_red(options.cas(ca_ix));
    tm_prob = p / sum(p(:));
    f(ca_ix,i) =  sum(tm_prob * (0:69)');   
  end
  nf(:,i) = f(:,i) / max(f(:,i));
end
figure(61); clf; hold on;
xsc = options.cas;
semilogx(xsc,nf(:,1),'-','Color',[0.1 0.1 0.1],'LineWidth',3);
semilogx(xsc,nf(:,2),'--','Color',[0.8 0.2 0.2],'LineWidth',2);
semilogx(xsc,nf(:,3),'-.','Color',[0.2 0.2 0.8],'LineWidth',2);
set(gca,'XScale','log'); xlim([1e-1 1e1]);

xlabel('Calcium concentration (\mu{}M)'); 
ylabel('Normalized force'); 
legend('Default simulation', '(1) Crossbridges blocked by blebbistatin', '(2) Crossbridges sensitized by dATP','Location','NorthOutside');
set(gca,'ytick',[0 0.25 0.5 0.75 1]); set(gca,'xtick',[1e-2 1e-1 1e-0 1e1])
myfig2(figure(61), './Fig6a.png',[3.25 2.5]*2,225);


 % rigor plot
 KDMS = 2.^(-8:0.01:8);
 options.cas = [0 1e6*10^-4.5];
 rnf=[];rf=[];
 
 for i=1:length(KDMS)
   KDM = KDMS(i);
   p_red = @(ca) SE .* (1 ./ options.KDA.^(options.n-nC_red)) .* (1 + ca/options.KDC).^(options.n-nC_red) .* (1 + ca/options.KDC + ca/(options.KDC*options.KDI)).^nC_red .* (1 ./ KDM.^nxb_red);
   for ca_ix = 1:length(options.cas)
     p = p_red(options.cas(ca_ix));
     tm_prob = p / sum(p(:));
     rf(ca_ix,i) =  sum(tm_prob * (0:69)');   
   end
 end   
 for i = 1:length(options.cas)
   rnf(i,:) = rf(i,:) / max(rf(i,:));
 end

 figure(62); clf; hold on;
 semilogx(KDMS,rnf(1,:),'--','Color',[0.2 0.2 0.8],'LineWidth',2);
 semilogx(KDMS,rnf(2,:),'-','Color',[0.8 0.2 0.2],'LineWidth',2);
 semilogx([options.KDM options.KDM],[0 1],':','Color',[0.1 0.1 0.1],'LineWidth',1.5);
 set(gca,'XScale','log');    
 xlabel('Myosin dissociation constant K_{DM}'); ylabel('Normalized force');
 legend('[Ca^{2+}] = 0 \mu{}M','[Ca^{2+}] = 31.6 \mu{}M (pCa 4.5)','Model default K_{DM} = 2','Location','NorthOutside');
 xlim([1e-2 1e2]);
 set(gca,'ytick',[0 0.25 0.5 0.75 1]); set(gca,'xtick',[1e-2 1e-1 1e-0 1e1 1e2])
   
 myfig2(figure(62), './Fig6b.png',[3.25 2.5750]*2,225);
 