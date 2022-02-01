# Chaotic Routes From Chains To Compact Structures
In this repository, I describe my most significant PhD result. 

## Background
In my Phd, I studied short chains of paramagnetic particles submerged in a fluid and subjected to an external rotating magnetic field. These particle systems have many proven and potential applications in microfluidics, optics and bio industry. I discovered in my simulations that above a critical rotational frequency of the field, the chain of the particles breaks up and particles undergo an episode of chaotic motion which is temporary. The long term response of the particles is forming ordered, rotating and stable structures. 

![](Figures/simulation.tif)

Then, we ran experiments and validated the simulation outcome: 

![](Figures/experiment.png)

## Tangent Dynamics Analyis 
In order to systematically capture the transition time from chaos to order I had to run tangent dynamics analysis. I do not go into the details of the analysis, however, I am sharing the code that re-generates the figure 4b of the paper which is shown below. A positive slope means chaos and a zero slope means order. The sharp change in slope is representative of the sharp transition from chaos to order. 

![](Figures/fig4b.png)

If you run the MATLAB code you will get the following figure which is the equivalent of figure 4b for a slighly different initial condition. 

![](Figures/fig4b_regenerated.png)

I have also attached a movie from real experiments that capture this sharp transition from chaos to order: experiment_chaos_order_transition.mp4

## Code
Please find the MATLAB code for the tangent dynamics analysis in the code folder. 


