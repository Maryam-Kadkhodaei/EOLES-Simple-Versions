$OnText
  A simple version of the EOLES model has been developed by me to undrestand the whole sysmtem better
$Offtext


sets
    h                                        /0*8759/
   first(h)        'first hour'
   last(h)         'last hour'
   tec             'Technology'             / biogas, pv, battery/
   gen(tec)             'energy generation tec'  /biogas, pv/
   comb            'comb technology'        /biogas/
   vre(tec)        'Renewable energy tec'   /pv/
   str(tec)        'Storage technology'     /battery/
;

first(h) = ord(h)=1;
last(h) = ord(h)=card(h);

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
$include  inputs/annuities - BiogasPVstr.csv
$offdelim
/;

parameter capa_max(tec) 'maximum capacities of the technologies in GW'
/
$ondelim
$include  inputs/max_capas - BiogasPVstr.csv
$offdelim
/ ;

parameter capex_en(str) 'annualized energy capex cost of storage technologies in M€/GWh/year'
/
$ondelim
$include  inputs/str_annuities - BiogasPVstr.csv
$offdelim
/ ;

parameter fOM(tec) 'annualized fixed operation and maintenance costs M€/GW/year'
/
$ondelim
$include  inputs/fO&M - BiogasPVstr.csv
$offdelim
/ ;
Parameter vOM(tec) 'Variable operation and maintenance costs in M€/GWh'
/
$ondelim
$include  inputs/vO&M - BiogasPVstr.csv
$offdelim
/ ;


*Parametres---------------------------------------------------------------------
parameter eta_in(str) 'charging efifciency of storage technologies' /battery 0.9/;
parameter eta_out(str) 'discharging efficiency of storage technolgoies' /battery 0.95/;
parameter s_capex(str) 'charging related annuity of storage in M€/GW/year' / battery 0/;
parameter s_opex(str)    'charging related fOM of storage in M€/GW/year'   / battery 0/;

* Model----------------------------------------------------------------------------
* Variable-------------------------------------------------------------------------
Variables
               GENE(tec,h)        'Generation of each technology at each hour'
               CAPA(tec)          'Capacity of each technology'
               STORAGE(str,h)  'hourly electricity input of battery storage GW'
               S(str)          'charging power capacity of each storage technology'
               CAPACITY(str)   'energy volume of storage technologies in GWh'
               STORED(str,h)   'energy stored in each storage technology in GWh'
               COST               'Total Cost'

positive variables GENE(tec,h),CAPA(tec),STORAGE(str,h), S(str),STORED(str,h),CAPACITY(str);

equations
             adequecy                  'Energy generation should satisfy the Demand'
             obj                       'The total cost of the energy'
             generation_limit          'Generation limit'
             gene_vre                  'VRE Sources energy generation'
             storage_const             'Energy stored is the same in the forst and last time step'
             storage_capa1             'the capacity with hourly charging relationship of storage'
             storage_capa2             'storage power limit'
             stored_cap                'maximum energy that is stored in storage units'
             stored_energy_battery     'Energy Stored in Battery';

gene_vre(vre,h)..                     GENE(vre,h)             =e=     CAPA(vre)*load_factor(vre,h);
storage_const(str, first, last)..     STORED(str,first)       =e=     STORED(str,last)+STORAGE(str,last)*eta_in(str)-GENE(str,last)/eta_out(str);
stored_energy_battery(h+1,h,str)..    STORED(str,h+1)         =e=     STORED(str,h) + STORAGE(str,h)*eta_in(str) -GENE(str,h)/eta_out(str);
stored_cap(str,h)..                   STORED(str,h)           =l=     CAPACITY(str);
storage_capa1(str,h)..                S(str)                  =g=     STORAGE(str,h);
storage_capa2(str)..                  S(str)                  =l=     CAPA(str);
generation_limit(tec,h)..             CAPA(tec)               =g=     GENE(tec,h);
adequecy(h)..                         GENE('pv',h)+ GENE('battery',h)+ GENE  ('biogas',h)  =g=     demand(h)+ STORAGE('battery',h) ;
obj..                                 COST                    =e=    (sum(tec,CAPA(tec)*capex(tec) ) + sum(str,CAPACITY(str)*capex_en(str)) + sum(str,S(str)*(s_capex(str)+s_opex(str)))+sum(tec,CAPA(tec)*fOM(tec)) + sum((tec,h),GENE(tec,h)*vOM(tec)))/1000;




*                                Model options
*-------------------------------------------------------------------------------
model RES_FR /all/;
*-------------------------------------------------------------------------------
*                                Solve statement
*-------------------------------------------------------------------------------
$If exist RES_FR_p.gdx execute_loadpoint 'RES_FR_p';
Solve RES_FR using lp minimizing COST;


CAPA.up('pv') = capa_max('pv');

parameter sumdemand      'the whole demand per year in TWh';
sumdemand =  sum(h,demand(h))/1000;

parameter gene_tec(tec) 'Overall yearly energy generated by the technology in TWh';
gene_tec(tec) = sum(h,GENE.l(tec,h))/1000;

Parameter lcoe(gen);
lcoe(gen) = (CAPA.l(gen)*(fOM(gen)+capex(gen))+ gene_tec(gen)*vOM(gen)*1000)/gene_tec(gen);

parameter lcos(str);
lcos(str) = (CAPA.l(str)*(fOM(str)+capex(str))+ gene_tec(str)*vOM(str) + S.l(str)*(s_capex(str)+s_opex(str))+ CAPACITY.l(str)*capex_en(str))/gene_tec(str);


parameter gene_vre_mult(vre,h) 'vre generation'    ;
gene_vre_mult(vre,h)= 1000000*gene_vre.m(vre,h);

*parameter storage_const_mult(str) 'storage constant'    ;
*storage_const_mult(str)= 1000000*storage_const.m(str,h);

parameter stored_energy_battery_mult(h,str) 'SoC'    ;
stored_energy_battery_mult(h,str)= 1000000*stored_energy_battery.m(h+1,h,str);


parameter str_power_lim_low(str,h) 'storage power limit lower'    ;
str_power_lim_low(str,h) = 1000000*storage_capa1.m(str,h);

parameter str_power_lim_high(str) 'storage power limit higher'    ;
str_power_lim_high(str) = 1000000*storage_capa2.m(str);

parameter store_cap_mult(str,h) 'storage '    ;
store_cap_mult(str,h) = 1000000*stored_cap.m(str,h);

parameter gen_lim_mult(tec,h) 'generation limit'    ;
gen_lim_mult(tec,h) = 1000000*generation_limit.m(tec,h);


parameter spot_price(h) 'marginal cost'    ;
spot_price(h) = 1000000*adequecy.m(h);

display lcos
display lcoe
file hourly_generation /'BiogasPVstrmult2_Output.csv' / ;
*the .csv file
put hourly_generation;
hourly_generation.pc=5;
put 'hour'; loop(tec, put tec.tl;); put 'Charging Power',put 'SoC' ,put 'demand' , 'price','SoC mult','gene_vre_mult', 'mul_pv','mul_biogas ','mul_battery','Store_Capacity','str_power_lim_low','str_power_lim_high'/ ;
loop (h,
put h.tl; loop(tec,  put GENE.l(tec,h);) put STORAGE.l('battery',h),put STORED.l('battery',h) ,put demand(h), spot_price(h),abs(stored_energy_battery_mult(h,'battery')),abs(gene_vre_mult('pv',h)),abs(gen_lim_mult('pv',h)),abs(gen_lim_mult('biogas',h)),abs(gen_lim_mult('battery',h)),abs( store_cap_mult('battery',h)),abs(str_power_lim_low('battery',h)), abs(str_power_lim_high('battery'))/
;);
