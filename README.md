# Team Biowhackathon Repository

Our team proposes a novel method of producing highly polarised magnets in microgravity, using a combination of biological agents as molecular scaffolding, and soundwaves to control the compartmentalisation and mixing of our reagents in a single vessel.

This repository demos a simple simulation

## Rationale



## Computational Fluid Dynamics (CMD) Simulation

This uses the Navier-Stokes method to simulate fluid dynamics in microgravity. This method treats the solvent as a continous medium, simulating motion as velocity fields rather than vibration of particles. Compared to simulating solvent as a particle, this method is computationally inexpensive, allowing for easy scaling up for clusters or scaling down for laptops (and anything in between). Tradeoffs of simulation accuracy for compute power is mainly decided by the number of nodes present and the length of timesteps. While there is not any interface for changing these, the variables 

