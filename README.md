# Modelling COVID-19 in the North American region with a metapopulation network and Kalman filter
medRxiv preprint: https://doi.org/10.1101/2024.06.05.24308495
_Metapopulation SEI<sup>r</sup>I<sup>u</sup>R model for North American region (CA, US, MX) at the state level. <br/>
Implemented with EAKF, daytime and nighttime transmission, daily adjusted mobility rate, integrated stochastically_
<p float="left">
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
