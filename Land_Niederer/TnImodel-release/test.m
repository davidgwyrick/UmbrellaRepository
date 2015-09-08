%% Script to use Land/Niederer functions to create Ca Transient from 2015 paper
% 9/8/15
%%
options = dynamic_xb_prep('human');
TopDir=pwd;
FigureDir=['J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\Land_Niederer\TnImodel-release\figures\'];
cd(FigureDir);
ca = make_ca(0.14,.3431,48.2,175.9,343.1);
cd(TopDir);
[t,y] = ode15s( @(t,y) make_dy(t,y,ca,options), 0:1000, options.y0);
dydt = make_dy(t,y,ca,options);
