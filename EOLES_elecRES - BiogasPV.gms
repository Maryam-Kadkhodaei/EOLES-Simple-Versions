$OnText
  A simple version of the EOLES model has been developed by me to undrestand the whole sysmtem better
$Offtext


sets
h                                        /0*8759/
tec             'Technology'             / biogas, pv/
vre(tec)           'Renewable energy tec'   /pv/

* Input data--------------------------------------------------------------------
$Offlisting
parameter load_factor(vre,h) 'Production profiles of VRE'
/
$ondelim
$include  inputs/vre_pv_2007.csv
$offdelim
/;

parameter demand(h) 'demand profile in each hour in GW'
/
$ondelim
$include inputs/demand2050_ademe.csv
$offdelim
/;

parameter capex(tec) 'annualized power capex cost in M€/GW/year'
/
$ondelim
$include  inputs/annuities - BiogasPV.csv
$offdelim
/;

parameter fOM(tec) 'annualized fixed operation and maintenance costs M€/GW/year'
/
$ondelim
$include  inputs/fO&M - BiogasPV.csv
$offdelim
/ ;
Parameter vOM(tec) 'Variable operation and maintenance costs in M€/GWh'
/
$ondelim
$include  inputs/vO&M - BiogasPV.csv
$offdelim
/ ;


* Model----------------------------------------------------------------------------
* Variable-------------------------------------------------------------------------
Variables
               GENE(tec,h)        'Generation of each technology at each hour'
               CAPA(tec)          'Capacity of each technology'
               COST               'Total Cost'

positive variables GENE(tec,h),CAPA(tec);

equations
             adequecy            'Energy generation should satisfy the Demand'
             obj                 'The total cost of the energy'
             generation_limit    'Generation limit'
             gene_vre            'VRE Sources energy generation' ;

gene_vre(vre,h)..                 GENE(vre,h)  =e=  CAPA(vre)*load_factor(vre,h);
generation_limit(tec,h)..         CAPA(tec) =g= GENE(tec,h);
adequecy(h)..                     sum(tec,GENE(tec,h)) =g= demand(h)  ;
obj..                             COST  =e= (sum(tec,CAPA(tec)*capex(tec) )+sum(tec,CAPA(tec)*fOM(tec)) + sum((tec,h),GENE(tec,h)*vOM(tec)))/1000;




*                                Model options
*-------------------------------------------------------------------------------
model RES_FR /all/;
*-------------------------------------------------------------------------------
*                                Solve statement
*-------------------------------------------------------------------------------
$If exist RES_FR_p.gdx execute_loadpoint 'RES_FR_p';
Solve RES_FR using lp minimizing COST;

parameter gene_vre_mult(vre,h) 'generation limit'    ;
gene_vre_mult(vre,h) = 1000000*gene_vre.m(vre,h);

parameter gen_lim_mult(tec,h) 'generation limit'    ;
gen_lim_mult(tec,h) = 1000000*generation_limit.m(tec,h);

parameter spot_price(h) 'marginal cost'    ;
spot_price(h) = 1000000*adequecy.m(h);


file hourly_generation /'BiogasPVmult2_Output.csv' / ;
*the .csv file
put hourly_generation;
hourly_generation.pc=5;
put 'hour'; loop(tec, put tec.tl;); put 'demand' , 'price', 'gen_lim_multpv','gen_lim_multbiiogas','gene_vre_mult'/ ;
loop (h,
put h.tl; loop(tec,  put GENE.l(tec,h);) put demand(h), spot_price(h),gen_lim_mult('pv',h),gen_lim_mult('biogas',h),gene_vre_mult('pv',h) /
;);
