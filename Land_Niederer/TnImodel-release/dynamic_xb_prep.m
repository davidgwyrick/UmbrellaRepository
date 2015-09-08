function options = dynamic_xb_prep(varargin)
 options = set_options(varargin{:}); %#ok<NOPRT>

 n=options.n;

 sefile = sprintf('SE%d_%.1f.mat',n,options.gamma);
 try
   load(sefile,'SE');
 catch err
   err
   error('SE file "%s" for n=%d, gamma=%.2f not found. Only n=26/gamma=2 is provided by default, please contact me for other files.',sefile,n,options.gamma);
 end
 
 num_xb = size(SE,2)-1;
 options.SE = SE;
 
 SEthr = options.SEthr; % limits ratios and rate constants a bit and reduce number of ODE states
 SE_small = SE./repmat(sum(SE,2),1,size(SE,2)) < SEthr;
 for i=1:size(SE,1)
    SE(i, find(SE_small(i,:)==0,1,'last')+1:end ) = 0;
 end
 
 [rst_i,xb_i] = find(SE ~= 0);
 SEix = SE; SEix( sub2ind(size(SE),rst_i,xb_i) ) = 1:length(rst_i); % number states
 
 totst = length(rst_i);
 %fprintf('Number of xb/tm states: %d\n',totst);
 
 kbc        = zeros(n+1,num_xb+1);
 kcb_no_tni = zeros(n+1,num_xb+1);
 xb_on      = zeros(n+1,num_xb+1);
 xb_off     = zeros(n+1,num_xb+1);
  
 clear j
 for i=0:n
  for j=0:num_xb
             SEi    = SE(i+1,j+1);            % prevent off-by-one
    if i~=n; SEplus = SE(i+2,j+1); kbc(i+1,j+1)       =  options.k_Aoff * (n-i) * (  (SEplus / (n-i)  ) / (SEi  / (i+1))  )^options.r     ; end
    if i~=0; SEmin  = SE(i,j+1);   kcb_no_tni(i+1,j+1)=  options.k_Aon  *  i    * (  (SEi    / (n-i+1)) / (SEmin/   i  )  )^-(1-options.r); end

    if j~=num_xb; SEplus = SE(i+1,j+2); xb_on(i+1,j+1)  = options.k_Xon  * (num_xb-j) * ( (SEplus/(num_xb-j)) / (SEi /(j+1)  ) )^options.q; end
    if j~=0;      SEmin  = SE(i+1,j);   xb_off(i+1,j+1) = options.k_Xoff * j          * ( (SEi / (num_xb-j+1) / (SEmin/ j )  ) )^-(1-options.q); end
  end
 end
 
 % fix any overflow, which can happen if SEthr=0
 xb_off(isinf(xb_off)) = 1e9; 
 xb_on(isnan(xb_on))            = 0;  options.xb_on     = xb_on;
 xb_off(isnan(xb_off)) = 0;           options.xb_off    = xb_off;
                                      
 kbc(isinf(kbc)) = 1e9;
 kcb_no_tni(isinf(kcb_no_tni)) = 1e9;
 kbc(isnan(kbc))                = 0;  options.kbc       = kbc;
 kcb_no_tni(isnan(kcb_no_tni))  = 0;  options.kcb_no_tni= kcb_no_tni;
 
 options.num_xb = num_xb;
 options.num_tmxb  = length(rst_i);
 options.state_tmi = rst_i;
 options.state_xbj = xb_i;
 options.state_num = SEix;
 options.frac_c    = (rst_i-1)/n;
 
 % initial condition
 y0 = zeros(options.num_tmxb+3,1); 
 y0(options.state_num(1,1)) = 1;
 options.y0=y0;
 
 
