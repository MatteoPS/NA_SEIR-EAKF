# Modelling COVID-19 in the North American region with a metapopulation network and Kalman filter<br/><br/>

_Metapopulation SEI<sup>r</sup>I<sup>u</sup>R model for North American region (CA, US, MX) at the state level. <br/>
Implemented with EAKF, daytime and nighttime transmission, daily adjusted mobility rate, integrated stochastically_
<p float="left">
medRxiv preprint: https://doi.org/10.1101/2024.06.05.24308495<br/><br/>
<br/><br/>

**Estimated parameters:**<br/>
**_α_** = ascertainment rate<br/>
**_β_** = transmission rate<br/>
<br/>
**Fixed parameters:**<br/>
**_Z_** = latency period<br/>
**_D_** = infectious period<br/>
**_µ_** = relative transmissibility<br/>
**_θ_** = mobility factor<br/><br/>

<br/><br/>
### How to run:
```Matlab
MODEL_RUN("run_name")
```
this funtion will run the model with 300 ensemble members to estimate **_α_** (ascertainment rate) and **_β_** (transmission rate) for the 96 location of the North American region.<br/><br/>
the function automatically call the plotting function to plot state variables and parameters for 9 selected locations.<br/><br/>
### Commuting matirix:
The daily commuting matrix was retrieved and derived form four datasest:<br/>
1) [Canadian 2016 census (Statistics Canada) Commuting Flow from Geography of Residence to Geography of Work](https://www12.statcan.gc.ca/census-recensement/2016/dp-pd/dt-td/Rp-eng.cfm?TABID=4&LANG=E&A=R&APATH=3&DETAIL=0&DIM=0&FL=A&FREE=0&GC=0&GL=-1&GID=1354564&GK=0&GRP=1&O=D&PID=111333&PRID=10&PTYPE=109445&S=0&SHOWALL=0&SUB=0&Temporal=2017&THEME=125&VID=0&VNAMEE=&VNAMEF=%20(2017)&D1=0&D2=0&D3=0&D4=0&D5=0&D6=0)<br/>
2) [Canada Frontier Counts (Statistics Canada): Number of vehicles travelling between Canada and the United States](https://www150.statcan.gc.ca/t1/tbl1/en/tv.action?pid=2410000201)<br/>
3) [5-Year American Community Survey (ACS) Commuting Flows (United States Census Bureau](https://www.census.gov/data/tables/2015/demo/metro-micro/commuting-flows-2015.html)<br/>
4) [Mexican Intercensal Survey 2015 (National Institute of Statistics and Geography, INEGI)](https://en.www.inegi.org.mx/programas/intercensal/2015/#Microdatos)<br/>
<br/>
The resulting commuting matrix contains the number of people that commute daily to work in another location
It's available here:<br/>
[Create_nl_part_Cave/final_commuting_matrix_Oct2023.csv](Create_nl_part_Cave/final_commuting_matrix_Oct2023.csv)

