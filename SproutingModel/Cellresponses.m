close all;
clear;


GoodFits = [1:18];

for Num_Fit = 1:size(GoodFits,2) 
    
load('rg_V.mat','rg_V');
load('rm_V.mat','rm_V');
load('P_V.mat','P_V');


Heiss_VEGF = [0	2	4	8	16	32	64
403.9270687	641.7952314	747.2650771	985.1332398	1052.454418	1045.7223	1115.287518
80.78541374	195.2314165	201.9635344	80.78541374	175.0350631	199.7194951	293.9691445];

InitCellNum_Heiss = 500;

avgHUVEC_diameter = 2*((1e3*3/(4*pi))^(1/3));
r_cell = (1e3/(pi*4/3))^(1/3);
R_sphere = (InitCellNum_Heiss*(r_cell^2)/4)^(1/2);
Numcell_R = floor(2*pi*R_sphere/(2*r_cell));

Ncell = ones(InitCellNum_Heiss, 1);
        


[a,b]= size(Ncell);
VEGFdose = Heiss_VEGF(1,:);
        

 for RepeatTimes = 1:10 
  for m = 1:size(VEGFdose,2)
    T_V = 24; % hours
    tstep = 1; % hours
    
V_avg(m,RepeatTimes) = rm_V(m,GoodFits(Num_Fit)); %um/h
DT_avg(m,RepeatTimes) = 1/rg_V(m,GoodFits(Num_Fit)); %h

Probcutoff_avg(m,RepeatTimes) = 1-P_V(m,GoodFits(Num_Fit));

sproutingdetail{m,RepeatTimes}(:,1) = 1:1:InitCellNum_Heiss; % the number of cell
sproutingdetail{m,RepeatTimes}(:,2) = 1:1:InitCellNum_Heiss; % the number of the mother cell


SD_DT = 12;
SD_V = 12;
SD_Prob = 3/(1-Probcutoff_avg(m,RepeatTimes)); % mean +3*SD < 1 -> SD < (1-mean)/3

for i = 1:a
DT(i,1) = findParamValue(DT_avg(m,RepeatTimes),SD_DT);
V(i,1) = findParamValue(V_avg(m,RepeatTimes),SD_V);
Probcutoff(i,1) = findParamValue(Probcutoff_avg(m,RepeatTimes),SD_Prob);

if (Probcutoff(i,1) >= 1) && (Probcutoff(i,1) <= 0)
    Probcutoff(i,1) = findParamValue(Probcutoff_avg,SD_Prob);
else
    Probcutoff(i,1) = Probcutoff(i,1);
end

sproutingdetail{m,RepeatTimes}(i,3) = DT(i,1); % cell doubling time
sproutingdetail{m,RepeatTimes}(i,4) = V(i,1); % sprout growth rate
sproutingdetail{m,RepeatTimes}(i,5) = Probcutoff(i,1); % probability cutoff of forming a sprout
end

sproutingdetail{m,RepeatTimes}(:,6) = 0; % 1 - tip cell, 0 - not a tip cell
sproutingdetail{m,RepeatTimes}(:,7) = 0; % time to become a tip cell
sproutingdetail{m,RepeatTimes}(:,8) = 0; % sprout length by the end of the simulation
sproutingdetail{m,RepeatTimes}(:,9) = 0; % the time that a new cell is generated
sproutingdetail{m,RepeatTimes}(:,10) = 0; % maximum number of tip cell can be formed at present



for t = 1:tstep:T_V
        for i = 1:a
            
            if sproutingdetail{m,RepeatTimes}(i,6) == 0
                
                P = rand;
                
                if P >= sproutingdetail{m,RepeatTimes}(i,5)
                    sproutingdetail{m,RepeatTimes}(i,6) = 1;
                    sproutingdetail{m,RepeatTimes}(i,7) = t;
                    sproutingdetail{m,RepeatTimes}(i,8) = (T_V-t)*sproutingdetail{m,RepeatTimes}(i,4);
                  sproutingdetail{m,RepeatTimes}(i,10) = 0; 
                  if t == T_V
                        sproutingdetail{m,RepeatTimes}(i,6) = 0;
                    end
                end
                
                if (P < sproutingdetail{m,RepeatTimes}(i,5)) && (sproutingdetail{m,RepeatTimes}(i,6) == 0)
                    Pmax = 10e-4;
                    r_cell = (1e3/(pi*4/3))^(1/3);
                    sproutingdetail{m,RepeatTimes}(i,10) = r_cell*2*Pmax*(t-sproutingdetail{m,RepeatTimes}(i,9));                    
                    
                    if (floor(t/sproutingdetail{m,RepeatTimes}(i,3)) - floor((t-1)/sproutingdetail{m,RepeatTimes}(i,3)) > 0)
                       sproutingdetail{m,RepeatTimes}(end+1,1) =  sproutingdetail{m,RepeatTimes}(end,1)+1;
                       sproutingdetail{m,RepeatTimes}(end,2) =  i;
                       sproutingdetail{m,RepeatTimes}(end,3) =  sproutingdetail{m,RepeatTimes}(i,3);
                       sproutingdetail{m,RepeatTimes}(end,4) =  sproutingdetail{m,RepeatTimes}(i,4);
                       sproutingdetail{m,RepeatTimes}(end,5) =  sproutingdetail{m,RepeatTimes}(i,5);
                       sproutingdetail{m,RepeatTimes}(end,6) =  sproutingdetail{m,RepeatTimes}(i,6);
                       sproutingdetail{m,RepeatTimes}(end,7) =  sproutingdetail{m,RepeatTimes}(i,7);
                       sproutingdetail{m,RepeatTimes}(end,8) =  sproutingdetail{m,RepeatTimes}(i,8);
                       sproutingdetail{m,RepeatTimes}(end,9) =  t;
                       sproutingdetail{m,RepeatTimes}(end,10) =  r_cell*2*Pmax*(t-sproutingdetail{m,RepeatTimes}(i,9));
                    end
                end
            end
            
                       
            
 N_tot(t) = size(sproutingdetail{m,RepeatTimes},1);
 N_tip(t) = sum(sproutingdetail{m,RepeatTimes}(:,6));
 
Tipcell_max(t) = floor(sum(sproutingdetail{m,RepeatTimes}(:,10)));
 Pmax = 10e-4;
  if N_tip(t) > Tipcell_max(t)
      sproutingdetail{m,RepeatTimes}(i,6) = 0;
sproutingdetail{m,RepeatTimes}(i,7) = 0;
sproutingdetail{m,RepeatTimes}(i,8) = 0;
sproutingdetail{m,RepeatTimes}(i,10) = r_cell*2*Pmax*(t-sproutingdetail{m,RepeatTimes}(i,9));                    
                    
  end
  
   N_tot(t) = size(sproutingdetail{m,RepeatTimes},1);
   N_tip(t) = sum(sproutingdetail{m,RepeatTimes}(:,6));
       

 end 
    
    
  
    a= size(sproutingdetail{m,RepeatTimes},1);
    
    
end
 

TotalTipCell(m,RepeatTimes) = N_tip(T_V);
TotalLength(m,RepeatTimes) = sum(sproutingdetail{m,RepeatTimes}(:,8));
if TotalTipCell(m,RepeatTimes) == 0
   AvgLength(m,RepeatTimes) = 0;
else
AvgLength(m,RepeatTimes) = TotalLength(m,RepeatTimes)/TotalTipCell(m,RepeatTimes);
end

[a,b]= size(Ncell);  


TotalTipCell_plane(m,RepeatTimes) = (pi/(2*((N_tot(T_V))^(1/2))))* TotalTipCell(m,RepeatTimes);
TotalLength_plane(m,RepeatTimes) = (pi/(2*((N_tot(T_V))^(1/2))))* TotalLength(m,RepeatTimes);

end 
  

  
 end


  if Num_Fit == 1
All_TotalTipCell = TotalTipCell_plane;
All_TotalLength = TotalLength_plane;
All_AvgLength = AvgLength;
  else
All_TotalTipCell = [All_TotalTipCell TotalTipCell_plane];
All_TotalLength = [All_TotalLength TotalLength_plane];
All_AvgLength = [All_AvgLength AvgLength];

  end
  

sproutingdetail =[];
end

All_TotalLength_VEGF = All_TotalLength;
All_TotalTipCell_VEGF = All_TotalTipCell;
All_AvgLength_VEGF = All_AvgLength;

save('All_TotalLength_VEGF.mat','All_TotalLength_VEGF');
save('All_TotalTipCell_VEGF.mat','All_TotalTipCell_VEGF');
save('All_AvgLength_VEGF.mat','All_AvgLength_VEGF');

figure(1)
subplot(1,3,1)
shadedErrorBar(VEGFdose,All_TotalLength',{@mean,@std},{'Color',[0.301960796117783 0.745098054409027 0.933333337306976],'LineWidth',2},2)

hXLabel = xlabel('VEGF (ng/ml)');
hYLabel = ylabel('Total Length (\mum)');
set([hXLabel, hYLabel]  , ...
'FontSize'   , 16          );

set(gca, ...
 'Box'         , 'off'     , ...   
 'TickDir'     , 'out'     , ...       % set the ticks to face outward
 'TickLength'  , [.015 .015] , ...     % set the length of the ticks
 'XMinorTick'  , 'on'      , ...       % show minor ticks
 'YMinorTick'  , 'on'      , ...
 'LineWidth'   , 2         ); 

subplot(1,3,2)
shadedErrorBar(VEGFdose,All_TotalTipCell',{@mean,@std},{'Color',[0.301960796117783 0.745098054409027 0.933333337306976],'LineWidth',2},2)

hXLabel = xlabel('VEGF (ng/ml)');
hYLabel = ylabel('Total Sprouts');
set([hXLabel, hYLabel]  , ...
'FontSize'   , 16          );

set(gca, ...
 'Box'         , 'off'     , ...   
 'TickDir'     , 'out'     , ...       % set the ticks to face outward
 'TickLength'  , [.015 .015] , ...     % set the length of the ticks
 'XMinorTick'  , 'on'      , ...       % show minor ticks
 'YMinorTick'  , 'on'      , ...
 'LineWidth'   , 2         ); 



subplot(1,3,3)
shadedErrorBar(VEGFdose,All_AvgLength',{@mean,@std},{'Color',[0.301960796117783 0.745098054409027 0.933333337306976],'LineWidth',2},2)

hXLabel = xlabel('VEGF (ng/ml)');
hYLabel = ylabel('Average length (\mum)');
set([hXLabel, hYLabel]  , ...
'FontSize'   , 16          );

set(gca, ...
 'Box'         , 'off'     , ...   
 'TickDir'     , 'out'     , ...       % set the ticks to face outward
 'TickLength'  , [.015 .015] , ...     % set the length of the ticks
 'XMinorTick'  , 'on'      , ...       % show minor ticks
 'YMinorTick'  , 'on'      , ...
 'LineWidth'   , 2         ); 

figure(2)

shadedErrorBar(VEGFdose,All_TotalLength',{@mean,@std},{'Color',[0.301960796117783 0.745098054409027 0.933333337306976],'LineWidth',2},2)
hold on
errorbar(Heiss_VEGF(1,:),Heiss_VEGF(2,:),Heiss_VEGF(3,:),  'MarkerSize',10,'Marker','diamond','LineStyle','none','linewidth', 2,'Color',[0, 0.4470, 0.7410])


hXLabel = xlabel('VEGF (ng/ml)');
hYLabel = ylabel('Total Length (\mum)');
set([hXLabel, hYLabel]  , ...
'FontSize'   , 16          );

set(gca, ...
 'Box'         , 'off'     , ...   
 'TickDir'     , 'out'     , ...       % set the ticks to face outward
 'TickLength'  , [.015 .015] , ...     % set the length of the ticks
 'XMinorTick'  , 'on'      , ...       % show minor ticks
 'YMinorTick'  , 'on'      , ...
 'LineWidth'   , 2         ); 