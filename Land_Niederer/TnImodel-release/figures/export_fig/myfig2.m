function myfig2(h,outfile,realsize,dpi,nowrite)
 if nargin<4; dpi=300; end
 if nargin<3; realsize=[3.25 2.5]; end

 ppi=get(0,             'ScreenPixelsPerInch');

 set(h,'position',[50 50 realsize* ppi])
 set(h,'color','w')

 ax=get(h,'CurrentAxes');
 set(ax,'fontname','Arial');

export_fig(outfile,'-png',sprintf('-r%d',dpi),'-painters','-nocrop',h);

