clear all
close all
type = 'weighted';

load('initvalue.mat','initvalue');
load('params.mat','params');

load('bestfits.mat','bestfits');

Heiss_VEGF = [0 2	4	8	16	32	64
403.927069 641.7952314	747.2650771	985.1332398	1052.454418	1045.7223	1115.287518];

    
for numfit = 1:size(bestfits,1)
POS = bestfits(numfit,:);
initvalue = initvalue(:,1);

%--- Assign the baseline parameter values 
params = params(:,1);

pos =POS;

V_max_pERK_m = pos(	1	);
V_max_pAkt_m = (48/65)*V_max_pERK_m;
V_max_pERK_prob = pos(	2	);
k = pos(	3	);
V_max_pAkt_prob = k*V_max_pERK_prob;
n_pERK_m = pos(	4	);
n_pERK_p = pos(	5	);
EC50_pERK_prob = [1513.70773408719];

basal_rm = pos(	6	);
basal_p = pos(	7	);

EC50_pERK=	[213144.644473928]; % median
EC50_pAkt=	[30143.8672732723]; % median
Vmax_p = [0.0581360421409674]; % median
Vmax_s = (28/50)*Vmax_p;
ng = [3.37571924458141]; % median
ns = ng;
basal_rg = [0.000618611431002887];% median


n_pAkt_m = n_pERK_m;
n_pAkt_p = n_pERK_p;

EC50_pAkt_m = [151245.978788799];
EC50_pAkt_prob = [46010.4360716057];
EC50_pERK_m = [124215.519121243];

V_max_pERK_g =	Vmax_p;
V_max_pAkt_g =	Vmax_s;
npERK_g =	ng;
npAkt_g =	ng;


InitCellNum = 500;
Ncell = ones(InitCellNum, 1);
avgHUVEC_diameter = 2*((1e3*3/(4*pi))^(1/3));

[a,b]= size(Ncell);
    % Set details
    options = odeset('RelTol',1e-9,'AbsTol',1e-12);

    Avogadro = 6.02214e23;

    Time = [0:60:120*60]; 
  
     % simulate the model
%% VEGF 
    
VEGFdose = Heiss_VEGF(1,:);
    % Specify initial values

  for m = 1:size(VEGFdose,2)
      
         if m ==1
           rg_V(m,numfit) = basal_rg;
           rm_V(m,numfit) = basal_rm;
           P_V(m,numfit) = basal_p* avgHUVEC_diameter;
         else
           
            VEGFconc = (7430000/50)*VEGFdose(m);
            initvalue(2,1) = 0; % FGF  
            initvalue(50,1) = VEGFconc;
            % simulate the model
            [~, predConcV] = ode15s(@coreFile_moleculardetail,Time,initvalue,options,params);
            pAkt_V(:,m) = sum(predConcV(:,[84,86,87,88,90,91,93]),2);
            maxpAkt_V(m,1) = max(pAkt_V(:,m));

            pERK_V(:,m) = sum(predConcV(:,[7,16,18,19,23,24,31,32]),2);
            maxpERK_V(m,1) = max(pERK_V(:,m));

            rg_pERK_V = V_max_pERK_g/(1+(EC50_pERK/maxpERK_V(m,1))^npERK_g);
            rg_pAkt_V = V_max_pAkt_g/(1+(EC50_pAkt/maxpAkt_V(m,1))^npAkt_g);          
            rg_V(m,numfit) = basal_rg + rg_pERK_V + rg_pAkt_V;

            rm_pERK_V = V_max_pERK_m/(1+(EC50_pERK_m/maxpERK_V(m,1))^n_pERK_m);
            rm_pAkt_V = V_max_pAkt_m/(1+(EC50_pAkt_m/maxpAkt_V(m,1))^n_pAkt_m);          
            rm_V(m,numfit) = basal_rm + rm_pERK_V + rm_pAkt_V;    

            P_pERK_V = V_max_pERK_prob/(1+(EC50_pERK_prob/maxpERK_V(m,1))^n_pERK_p);
            P_pAkt_V = V_max_pAkt_prob/(1+(EC50_pAkt_prob/maxpAkt_V(m,1))^n_pAkt_p);          
            P_V(m,numfit) = (basal_p + P_pERK_V + P_pAkt_V)* avgHUVEC_diameter;
    
         end
save('rg_V.mat','rg_V');
save('rm_V.mat','rm_V');
save('P_V.mat','P_V');
  end
end