# Coupled Wind Turbine Design and Layout Optimization
![Alt Text](/FLORISSE3D/src/FLORISSE3D/optAlturas/vid_f/small7b.gif)

In this study, we have discovered that wind farms can be significantly improved by designing wind turbines and the wind farms where they will operate at the same time! Specifically, there are many situations where having different wind turbine designs in the same farm results in higher energy production and lower cost of energy. Because of the difference in turbine height and rotor diameter, turbines in one group have minimal interactions with the turbines from the other group. 

The animation above shows the an example of the optimization process used in this research. This optimization used the wind rose shown on the bottom left, and assumed a constant velocity from each direction. As shown, the turbine positions (in the top left), heights, and rotor diameters (in the top middle) are all modified at the same time. What isn't seen in this example is that the tower diameter and thickness, turbine rating, and blade chord and twist distributions are also design variables. Also shown in the animation is the direction power production of the wind farm. Notice that in the optimal design there is significant power production from all all wind directions. The figures on the far right show the current cost of energy (on the top), and annual energy production (on the bottom). The optimal wind farm produces more energy, for a cheaper cost than the starting farm.

In this study we coupled wind farm power model, a wind farm cost model, a turbine tower structural model, and a surrogate of a complex rotor structural model in order to analyze wind farms. We considered wind speed variation with height from the atmospheric boundary layer, wind direction distributions from real wind farms, and a weibull wind speed distribution at every wind direction. We also calculated exact analytic gradients for every portion of the simulation to allow for good convergence of the many design variables during optimization.

Our results indicate that there can be significant benefits to building wind farms with different turbine designs. For wind farms with closely spaced with turbines, a wind farm with different turbine designs can have a cost of energy that is 7-10% lower than the same wind farm optimized with all the same turbines. 

Check out the full paper [here!](https://www.wind-energ-sci.net/4/99/2019/wes-4-99-2019.html)

Also take a look at another paper of ours that recently got published about wind farms with turbines that are different heights [here.](https://rdcu.be/blkdh)
