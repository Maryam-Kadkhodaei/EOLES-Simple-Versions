$OnText
  A simple version of the EOLES model has been developed by me to undrestand the whole sysmtem better
$Offtext


sets       h                                       /0*8759/
* tec             'technology'                    / biogas/

* Input data--------------------------------------------------------------------
parameter demand(h) 'demand profile in each hour in GW'
/
$ondelim
$include inputs/demand2050_ademe.csv
$offdelim
/;

parameter capex 'annualized power capex cost in M€/GW/year'
/
$ondelim
$include  inputs/annuities - biogas.csv
$offdelim
/;

parameter fOM 'annualized fixed operation and maintenance costs M€/GW/year'
/
$ondelim
$include  inputs/fO&M - biogas.csv
$offdelim
/ ;
Parameter vOM 'Variable operation and maintenance costs in M€/GWh'
/
$ondelim
$include  inputs/vO&M - biogas.csv
$offdelim
/ ;


* Model----------------------------------------------------------------------------
* Variable-------------------------------------------------------------------------
Variables
               GENE(h)        'Generation of each technology at each hour'
               CAPA          'Capacity of each technology'
               COST               'Total Cost'

positive variables GENE(h),CAPA;

equations
             adequecy        'Energy generation should satisfy the Demand'
             obj              'The total cost of the energy'
             generation_limit    'Generation limit';



obj..                         COST  =e= (CAPA*capex + CAPA*fOM + sum(h,GENE(h)*vOM))/1000;
adequecy(h)..                 GENE(h) =g= demand(h)  ;
generation_limit(h)..         CAPA =g= GENE(h);


*                                Model options
*-------------------------------------------------------------------------------
model RES_FR /all/;
*-------------------------------------------------------------------------------
*                                Solve statement
*-------------------------------------------------------------------------------
$If exist RES_FR_p.gdx execute_loadpoint 'RES_FR_p';
Solve RES_FR using lp minimizing COST;

parameter spot_price(h) 'marginal cost'    ;
spot_price(h) = 1000000*adequecy.m(h);

parameter gen_Limit_mult(h) 'gen_Limit_mult'    ;
gen_Limit_mult(h) = 1000000*generation_limit.m(h);

file hourly_generation /'Biogas_Output_mult.csv' / ;
*the .csv file

put hourly_generation;
hourly_generation.pc=5;
*put 'hour'; loop(h,); put 'demand' , 'price'/ ;
loop (h,
put h.tl;  put GENE.l(h), put demand(h), spot_price(h),gen_Limit_mult(h)/
;);

display CAPA.l