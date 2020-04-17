# CrystalAFC
MATLAB scripts to simulate anti-fouling control of plug-flow crystallization via heating and cooling cycle as described in https://doi.org/10.1016/j.ifacol.2015.08.180 and https://doi.org/10.1109/LLS.2017.2661981

Specifically, the scripts solve a set of partial differential equations that describes a crystallization and encrustation process in a plug flow:

![](Images/PFC-PBM_PDE.png)

Here, n is the crystal size distribution (CSD), u z is the slurry flow velocity, G is the crystal growth rate, B is the nucleation rate, nseed is the seed CSD, z is the reactor axis, L is the crystal size axis, and t is the time axis. A_{f}(t, z) = πR_{f}^{2}(t, z) is the flow area within the tube which changes with time and along the reactor length due to encrustation. A model for encrustation in a PFC inspired from fouling kinetics commonly found in heat exchangers is shown in the figure below. 

![](Images/PFC_domain.jpg)

The domains of the crystallization and encrustation are divided into three regions, namely the reactor wall (Ω_{W} : r ∈ [R_{f}, R_{0}]), the encrust (Ω_{E} : r ∈ [R_{i} , R_{f} ]) and the convection within the tube (Ω_{T} : r ∈ [0, R_{i} ]). Additionally, the encrustation dynamics can be summarized as follows:

![](Images/encrust_ODE.png)

Where m_{d} and m_{r} are the mass deposited and removed, respectively, δ is the encrust thickness, k_{E} is the thermal conductivity, χ is the thermal resistance, ρ_{E} is the encrust density, m is the encrust mass, k_{m} is the mass transfer coefficient of solute from the bulk solution to the encrust film, k_{R} is the adsorption rate of solute to encrust, C_{b} is the bulk solute concentration, C_{sat} is the saturation concentration within the boundary or film layer, w is the bulk fluid velocity. α is the linear expansion coefficient, ∆T is the temperature difference between the reactor wall and the encrust surface, d_{p} is the encrust particle diameter, η is the film viscosity, and g is the gravitational acceleration. The adsorption rate is modeled as an Arrhenius-type expression where kR0 is the adsorption rate constant, ∆Ef is the activation energy, R is the ideal gas constant, and T_{f} is the film temperature.

Please look into the articles for more details, including the mass and energy balances. 
