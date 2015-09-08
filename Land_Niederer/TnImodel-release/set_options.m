% Prepares option structure and sets parameters
% options = set_options(varargin)
%  varargin given as a single string
%   mouse/rat/human/mouse_nosplit/rat_nosplit/human_nosplit: pre-parametrized models (nosplit implies q=0, otherwise q=0.5)
%  varargin given as ('option',value,'option',value) pairs for model parameters
%    gamma = 2
%    KDM   = 2
%    KDA   = 1e-3
%    KDI   = 4e-3
%    KDC   = 5.9 (microM)
%    r     = 0.5
%    q     = 0.5
%    k_Aoff= 10   ("k_A-") on-rate set using KDA
%    k_Ioff= 10   ("k_I-") on-rate set using KDI
%    k_Xoff= 0.5  ("k_M-") on-rate set using KDM
%    k_Con = 0.5  ("k_C+") off-rate set using KDC
%  (as well as some model running options. TODO: filter and doc/remove)
function options = set_options(varargin)
  if nargin==1
    if strcmp(varargin{1},'mouse')
      options = set_options(); return;
    elseif strcmp(varargin{1},'rat')
      options = set_options('k_Xoff',0.09,'KDI',0.0087); return;
    elseif strcmp(varargin{1},'human')
      options = set_options('n',26,'k_Xoff',0.016,'KDI',3.3e-3); return; % make_ca based
    elseif strcmp(varargin{1},'mouse_nosplit')
      options = set_options('k_Xoff',2,'KDI',4e-3,'q',1); return;
    elseif strcmp(varargin{1},'rat_nosplit')
      options = set_options('k_Xoff',0.3,'KDI',8.7e-3,'q',1); return;
    elseif strcmp(varargin{1},'human_nosplit')
      options = set_options('n',26,'k_Xoff',0.049,'KDI',3.3e-3,'q',1); return; % make_ca based
    end
    error('Expecting mouse/rat/human[_nosplit] or name-value pairs!');
  end

  if mod(nargin,2)~=0; error('Expecting name-value pairs!'); end
  
  % default is mouse with q=0.5
  options.gamma = 2;
  options.KDM   = 2;  
  options.KDA   = 1e-3;
  options.KDI   = 4e-3;
  options.KDC   = 5.9;

  % ----------
  options.n  = 26;  
  options.r  = 0.5;
  options.q  = 0.5;

  options.k_Con = 0.5;  
  options.k_Aoff= 10;
  options.k_Ioff= 10;
  options.k_Xoff= 0.5;

  % model running options
  options.tss = 100; % time to ss
  options.debug = 0;
  options.quick = 0; % less info/results
  options.use_st = 1;
  options.SEthr = 1e-6; % SE / sum(SE row) < this,  is considered 0 in dynamic_xb
  options.pca   = 4:0.1:8;
  
  for i=1:2:nargin
    assert(isfield(options,varargin{i}),sprintf('Option %s not recognized!',varargin{i}));
    options.(varargin{i}) = varargin{i+1};
  end
  
  options.k_Aon  = options.k_Aoff / options.KDA;
  options.k_Ion  = options.k_Ioff  / options.KDI;
  options.k_Xon  = options.k_Xoff / options.KDM;
  options.k_Coff = options.k_Con   * options.KDC; 
  
  options.cas =  10.^-options.pca(end:-1:1) * 1e6;

  
