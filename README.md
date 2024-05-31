# North-America-Metapop

_Metapopulation SEI<sup>r</sup>I<sup>u</sup>R model for North American region (CA, US, MX) at the state level. <br/>
Implemented with EAKF, reporting delay, daytime and nighttime transmission, daily adjusted mobility rate, integrated stochastically_
<p float="left">
<img src=https://github.com/ShamanLabDev/North-America-Metapop/assets/32901863/050194a0-446f-46e1-b7fa-4b0a6f0010e4 width="30%" />
<img src=https://github.com/ShamanLabDev/North-America-Metapop/assets/32901863/c13b94bb-9ef6-41fa-941a-ba1c1ffc3d3e width="55%" />
</p>
<br/><br/>

**Estimated parameters:**<br/>
**_α_** = reporting rate<br/>
**_β_** = transmission rate<br/>
<br/>
**Fixed parameters:**<br/>
**_Z_** = latency period<br/>
**_D_** = infectious period<br/>
**_µ_** = relative transmissibility<br/>
**_θ_** = mobility factor<br/><br/>
**parameters missing in this versions:**<br/>
**_L_** = loss of immunity<br/>
**_V_** = vaccination rate<br/>
<br/><br/>
### How to run (syntax not updated):
lunch a single run from matlab console
```Matlab
SIMULATIONS(num_ens, truth_name, run_nickname)
SIMULATIONS_glob_alphabeta(num_ens, truth_name, run_nickname)
%example run
SIMULATIONS_glob_alphabeta(150,"truth_NY_alpha02.mat","a02_glob","resample")
```
lunch multiple runs form bash using ``RUN_simulations.m``
```bash
/Applications/MATLAB_R2022a.app/bin/matlab -nodisplay -nosplash -nodesktop -r "run('RUN_simulations.m');exit;" | tail -n +11 
```

<br/><br/><br/><br/>
### Main functions:
- Estimating **_α_** and **_β_** as global parameters, updating the EAKF with the average of the adjustement from each state
  ```bash
  SIMULATIONS_glob_alphabeta.m
  ```
<br/><br/>
- Estimating **_α_** and **_β_** as local parameters
  ```bash
  SIMULATIONS.m
  ```
<br/><br/>
- Creates a synthetic truth and adds random noise to the observed variable. <br/>
  The output of this script (e.g. ``truth*.mat``) is required for any ``SIMULATIONS*.m`` to run 
  ```bash
  TRUTH.m
  ```
  
  <br/><br/>
- Plots state variables, parameters and EAKF variables for selected states. <br/>
  ```bash
  Plotting.m
  ```
  It is automaticcally called by ``SIMULATIONS*.m`` but it can be called using the output ``run_name.mat`` file saved at the end of  ``SIMULATIONS*.m``
  ```Matlab
  Plotting("run_name.mat")
  ```
