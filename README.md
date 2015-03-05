# Lightweight versus Efficiency under Cost Constraints.

This MATLAB codebase supports the paper: "APPLYING OPTIMAL CHOICES FOR REAL POWERTRAIN AND LIGHTWEIGHTING TECHNOLOGY OPTIONS TO PASSENGER VE-HICLES UNDER UNCERTAINTY" Accepted, Transport 2015 if you find it useful in any way, please cite the paper as such:

<TBC>

The files should run without any special MATLAB toolboxes.

The optimization routine is run by 'LWcase_knapcask' which is based on 

http://www.mathworks.com/matlabcentral/fileexchange/22783-0-1-knapsack/content/knapsack.m by Petter Strandmark

and is called by the scripts described below:

LWcase_no_sens – no sensitivity analysis performed

LWcase_knapsack_lifetime_sens1&2 - lifetime sensitivity analsis (revisions)

LWcase_knapsack_lw_sens- lightweighting technology cost sensitivity

LWcase_knapsack_eff_sens - efficiency technology cost sensitivity

LWcase_knapsack_fuel_sens - fuel cost sensitivity (forked from 4)

LWcase_knapsack_km_sens – driven km sensitivity (forked from 7)
