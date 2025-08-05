# Modeling COVID-19 in the North American Region with a Metapopulation Network and Kalman Filter

_Metapopulation SEI<sup>r</sup>I<sup>u</sup>R model for the North American region (CA, US, MX) at the state level._  
_Implemented with EAKF, incorporating daytime and nighttime transmission, daily adjusted mobility rates, and stochastic integration._  

**Publication** [https://doi.org/10.1016/j.epidem.2025.100818](https://doi.org/10.1016/j.epidem.2025.100818)

![Map of Commuting Flows](https://github.com/user-attachments/assets/aef9b42d-e2cb-4355-87fc-d3d9e670cbb0)  
![Metapopulation Structure](https://github.com/user-attachments/assets/2f386e76-9c25-42fe-9cb7-76be5f799259)  



### Estimated Parameters  
**_α_** = Ascertainment rate  
**_β_** = Transmission rate  

### Fixed Parameters  
**_Z_** = Latency period  
**_D_** = Infectious period  
**_µ_** = Relative transmissibility  
**_θ_** = Mobility factor  

---

### How to Run the Model  

```matlab
MODEL_RUN("run_name")
```
This function runs the model with 300 ensemble members to estimate **_α_** (ascertainment rate) and **_β_** (transmission rate) for the 96 locations in the North American region.  
It also automatically calls a plotting function to visualize state variables and parameters for 9 selected locations.  

---
### Commuting matrix:
The daily commuting matrix was retrieved and derived from four datasets:  

1) [Canadian 2016 census (Statistics Canada): Commuting Flow from Geography of Residence to Geography of Work](https://www12.statcan.gc.ca/census-recensement/2016/dp-pd/dt-td/Rp-eng.cfm?TABID=4&LANG=E&A=R&APATH=3&DETAIL=0&DIM=0&FL=A&FREE=0&GC=0&GL=-1&GID=1354564&GK=0&GRP=1&O=D&PID=111333&PRID=10&PTYPE=109445&S=0&SHOWALL=0&SUB=0&Temporal=2017&THEME=125&VID=0&VNAMEE=&VNAMEF=%20(2017)&D1=0&D2=0&D3=0&D4=0&D5=0&D6=0)  
2) [Canada Frontier Counts (Statistics Canada): Number of vehicles travelling between Canada and the United States](https://www150.statcan.gc.ca/t1/tbl1/en/tv.action?pid=2410000201)  
3) [5-Year American Community Survey (ACS) Commuting Flows (United States Census Bureau)](https://www.census.gov/data/tables/2015/demo/metro-micro/commuting-flows-2015.html)  
4) [Mexican Intercensal Survey 2015 (National Institute of Statistics and Geography, INEGI)](https://en.www.inegi.org.mx/programas/intercensal/2015/#Microdatos)  

The resulting commuting matrix contains the number of people that commute daily to work in another location.  

It's available as a CSV file here:  
[Create_nl_part_Cave/final_commuting_matrix_Oct2023.csv](Create_nl_part_Cave/final_commuting_matrix_Oct2023.csv)
