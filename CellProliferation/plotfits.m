close all;
clear;

load('besterrors.mat');
load('bestfits.mat');
load('best_r_p.mat');
load('best_r_s.mat');
load('best_r_g.mat');
load('best_r_p_V.mat');
load('best_r_s_V.mat');
load('best_r_g_V.mat');

r_p = best_r_p;
r_s = best_r_s;
r_g = best_r_g;
r_p_V = best_r_p_V;
r_s_V = best_r_s_V;
r_g_V = best_r_g_V;

GoodFits = 1:size(bestfits,1);
load('initvalue.mat','initvalue');
load('params.mat','params');

%  experimental data 
West_FGF = [0	0.039948403	0.040518737	0.041089071
0.01	0.017220343	0.017220876	0.017221409
0.03	0.015005104	0.016713972	0.018422841
0.1	0.025894211	0.052674305	0.079454399
0.3	0.064704613	0.099462447	0.134220282
1	0.118899096	0.136562914	0.154226732
3	0.148591089	0.167964842	0.187338594
10	0.207345044	0.237543941	0.267742838
30	0.135614668	0.184616451	0.233618234];

Jih_VEGF = [0	0.781914894	1.013297872	1.244680851
0.1	0.877659574	1.061170213	1.244680851
0.5	0.917553191	1.14893617	1.380319149
1	1.460106383	1.643617021	1.82712766
2	1.515957447	1.819148936	2.122340426
5	1.675531915	1.731382979	1.787234043
10	1.507978723	1.579787234	1.651595745];

Jih_FGF = [0	1.013700041	1.013803049	1.013906057
0.5	1.357643181	1.706118665	2.054594149
1	1.896992171	1.930881747	1.964771323
1.25	2.138545529	2.537958385	2.93737124
2.5	1.862124021	1.878914297	1.895704574
5	1.82303255	1.916512155	2.009991759
10	1.835290482	1.95421302	2.073135558];



FGFdose_data_West = West_FGF(:,1);
data_ref_F_West = West_FGF(6,3);

VEGFdose_data_Jih = Jih_VEGF(:,1);
data_ref_V_Jih = Jih_VEGF(4,3);

FGFdose_data_Jih = Jih_FGF(:,1);
data_ref_F_Jih = Jih_FGF(3,3);


for i = 1:size(FGFdose_data_West,1) 
    
    data_relative_F_West(i,1) = (West_FGF(i,2) - data_ref_F_West)/data_ref_F_West;
    data_relative_F_West(i,2) = (West_FGF(i,3) - data_ref_F_West)/data_ref_F_West;
    data_relative_F_West(i,3) = (West_FGF(i,4) - data_ref_F_West)/data_ref_F_West;
    west_neg(i) = data_relative_F_West(i,1) - data_relative_F_West(i,2);
    west_pos(i) = data_relative_F_West(i,3) - data_relative_F_West(i,2);
end

for i = 1:size(VEGFdose_data_Jih,1) 
    
    data_relative_V_Jih(i,1) = (Jih_VEGF(i,2) - data_ref_V_Jih)/data_ref_V_Jih;
    data_relative_V_Jih(i,2) = (Jih_VEGF(i,3) - data_ref_V_Jih)/data_ref_V_Jih;
    data_relative_V_Jih(i,3) = (Jih_VEGF(i,4) - data_ref_V_Jih)/data_ref_V_Jih;
    V_Jih_neg(i) = data_relative_V_Jih(i,1) - data_relative_V_Jih(i,2);
    V_Jih_pos(i) = data_relative_V_Jih(i,3) - data_relative_V_Jih(i,2);
end

for i = 1:size(FGFdose_data_Jih,1) 
    
    data_relative_F_Jih(i,1) = (Jih_FGF(i,2) - data_ref_F_Jih)/data_ref_F_Jih;
    data_relative_F_Jih(i,2) = (Jih_FGF(i,3) - data_ref_F_Jih)/data_ref_F_Jih;
    data_relative_F_Jih(i,3) = (Jih_FGF(i,4) - data_ref_F_Jih)/data_ref_F_Jih;
    F_Jih_neg(i) = data_relative_F_Jih(i,1) - data_relative_F_Jih(i,2);
    F_Jih_pos(i) = data_relative_F_Jih(i,3) - data_relative_F_Jih(i,2);
end
    
    


for numfit = 1:size(GoodFits,2)  
    FGFdose = [0:0.01:0.1 0.2 0.3 0.4 0.5 1:30];
    for i = 1:size(FGFdose,2) 
        Ncell_int_F = 1e4;
        T_F = 2*24;
        for j = 1:Ncell_int_F
            r_g_dis5e3_F (i,j) = findParamValue(r_g(i,GoodFits(numfit)),12);


            Nt_F(i,j) = (2^(floor(T_F*r_g_dis5e3_F (i,j))));

        end

        Ntot_F(i) = sum(Nt_F(i,:));
    end
    
    Ntot_ref_F = Ntot_F(16);


    for i = 1:size(FGFdose,2) 
        Ntot_relative_F(numfit,i) = (Ntot_F(i) - Ntot_ref_F)/Ntot_ref_F;

    end

VEGFdose = [0:0.1:1 2:1:10];


   for m = 1:size(VEGFdose,2)
        Ncell_int_V = 5000;
    T_V = 2*24;
    for j = 1:Ncell_int_V
        r_g_dis5e3_V (m,j) = findParamValue(r_g_V(m,GoodFits(numfit)),12);

            
        Nt_V(m,j) = (2^(floor(T_V*r_g_dis5e3_V (m,j))));

    end 
    
    Ntot_V(m) = sum(Nt_V(m,:));


    end
    
Ntot_ref_V = Ntot_V(11);


for i = 1:size(VEGFdose,2) 
    Ntot_relative_V(numfit,i) = (Ntot_V(i) - Ntot_ref_V)/Ntot_ref_V;

       
end
 

 
 
figure(1)
subplot(1,2,1)
plot(FGFdose, Ntot_relative_F(numfit,:),'linewidth', 2)
hold on

subplot(1,2,2)
plot(VEGFdose, Ntot_relative_V(numfit,:),'linewidth', 2)
hold on

end
figure(1)
subplot(1,2,1)
errorbar(FGFdose_data_West(1:6)', data_relative_F_West(1:6,2)',  west_pos(1:6),'MarkerSize',10,'Marker','o','LineStyle','none','linewidth', 2,'Color',[0.749019622802734 0.749019622802734 0])
hold on
errorbar(FGFdose_data_West(7:9)', data_relative_F_West(7:9,2)',  west_pos(7:9),'MarkerSize',10,'Marker','o','LineStyle','none','linewidth', 2,'Color',[0.87058824300766 0.490196079015732 0])
hold on
errorbar(FGFdose_data_Jih', data_relative_F_Jih(:,2)',  F_Jih_pos,'MarkerSize',10,'Marker','square','LineStyle','none','linewidth', 2,'Color',[0.87058824300766 0.490196079015732 0])
    hXLabel = xlabel('FGF2 (ng/ml)');
    hYLabel = ylabel('relative proliferation');
    set([hXLabel, hYLabel]  , ...
    'FontSize'   , 16          );

    set(gca, ...
     'Box'         , 'off'     , ...   
     'TickDir'     , 'out'     , ...       % set the ticks to face outward
     'TickLength'  , [.015 .015] , ...     % set the length of the ticks
     'XMinorTick'  , 'on'      , ...       % show minor ticks
     'YMinorTick'  , 'on'      , ...
     'LineWidth'   , 2         ); 
 set(gca, 'XScale', 'log')
 
subplot(1,2,2)
errorbar(VEGFdose_data_Jih(1:4)', data_relative_V_Jih(1:4,2)',V_Jih_pos(1:4),'MarkerSize',10,'Marker','square','LineStyle','none','linewidth', 2,'Color',[0.301960796117783 0.745098054409027 0.933333337306976])
hold on
errorbar(VEGFdose_data_Jih(5:7)', data_relative_V_Jih(5:7,2)',V_Jih_pos(5:7), 'MarkerSize',10,'Marker','square','LineStyle','none','linewidth', 2,'Color',[0.30588236451149 0.396078437566757 0.580392181873322])

hXLabel = xlabel('VEGF (ng/ml)');
hYLabel = ylabel('relative proliferation');
set([hXLabel, hYLabel]  , ...
'FontSize'   , 16          );

set(gca, ...
 'Box'         , 'off'     , ...   
 'TickDir'     , 'out'     , ...       % set the ticks to face outward
 'TickLength'  , [.015 .015] , ...     % set the length of the ticks
 'XMinorTick'  , 'on'      , ...       % show minor ticks
 'YMinorTick'  , 'on'      , ...
 'LineWidth'   , 2         ); 
 set(gca, 'XScale', 'log')
 
figure(2)
subplot(1,2,1)
shadedErrorBar(FGFdose(1,2:end),Ntot_relative_F(:,2:end),{@mean,@std},{'Color',[0.9290, 0.6940, 0.1250],'LineWidth',2},2)
hold on
errorbar(FGFdose_data_West(1:6)', data_relative_F_West(1:6,2)',  west_pos(1:6),'MarkerSize',10,'Marker','o','LineStyle','none','linewidth', 2,'Color',[0.749019622802734 0.749019622802734 0])
hold on
errorbar(FGFdose_data_West(7:9)', data_relative_F_West(7:9,2)',  west_pos(7:9),'MarkerSize',10,'Marker','o','LineStyle','none','linewidth', 2,'Color',[0.87058824300766 0.490196079015732 0])
hold on
errorbar(FGFdose_data_Jih', data_relative_F_Jih(:,2)',  F_Jih_pos,'MarkerSize',10,'Marker','square','LineStyle','none','linewidth', 2,'Color',[0.87058824300766 0.490196079015732 0])
    hXLabel = xlabel('FGF2 (ng/ml)');
    hYLabel = ylabel('relative proliferation');
    set([hXLabel, hYLabel]  , ...
    'FontSize'   , 16          );

    set(gca, ...
     'Box'         , 'off'     , ...   
     'TickDir'     , 'out'     , ...       % set the ticks to face outward
     'TickLength'  , [.015 .015] , ...     % set the length of the ticks
     'XMinorTick'  , 'on'      , ...       % show minor ticks
     'YMinorTick'  , 'on'      , ...
     'LineWidth'   , 2         ); 
 set(gca, 'XScale', 'log')
%  axis([0.01 30 -1 1]);
 
subplot(1,2,2)
shadedErrorBar(VEGFdose(1,2:end),Ntot_relative_V(:,2:end),{@mean,@std},{'Color',[0, 0.4470, 0.7410],'LineWidth',2},2)
hold on;
errorbar(VEGFdose_data_Jih(1:4)', data_relative_V_Jih(1:4,2)',V_Jih_pos(1:4),'MarkerSize',10,'Marker','square','LineStyle','none','linewidth', 2,'Color',[0.301960796117783 0.745098054409027 0.933333337306976])
hold on
errorbar(VEGFdose_data_Jih(5:7)', data_relative_V_Jih(5:7,2)',V_Jih_pos(5:7), 'MarkerSize',10,'Marker','square','LineStyle','none','linewidth', 2,'Color',[0.30588236451149 0.396078437566757 0.580392181873322])

hXLabel = xlabel('VEGF (ng/ml)');
hYLabel = ylabel('relative proliferation');
set([hXLabel, hYLabel]  , ...
'FontSize'   , 16          );

set(gca, ...
 'Box'         , 'off'     , ...   
 'TickDir'     , 'out'     , ...       % set the ticks to face outward
 'TickLength'  , [.015 .015] , ...     % set the length of the ticks
 'XMinorTick'  , 'on'      , ...       % show minor ticks
 'YMinorTick'  , 'on'      , ...
 'LineWidth'   , 2         ); 
 set(gca, 'XScale', 'log')

