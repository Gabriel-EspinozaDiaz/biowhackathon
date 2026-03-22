# Team Biowhackathon Repository

Our team proposes a novel method of producing highly polarised magnets in microgravity, using a combination of biological agents as molecular scaffolding, and soundwaves to control the compartmentalisation and mixing of our reagents in a single vessel. Full details are available in the report_metal_nucleating_protein_design documents (in both pdf and md formats). This quick guide goes over what the two demo simulations show, and how to run them.

## Computational Fluid Dynamics (CMD) Simulation

This uses the Navier-Stokes method to simulate fluid dynamics in microgravity. This method treats the solvent as a continous medium, simulating motion as velocity fields rather than vibration of particles. Compared to simulating solvent as a particle, this method is computationally inexpensive, allowing for easy scaling up for clusters or scaling down for laptops (and anything in between). Tradeoffs of simulation accuracy for compute power is mainly decided by the number of nodes present and the length of timesteps. While there is not any interface for changing these, the variables, upscaling could do with more nodes and shorter timesteps to get a simulation closer to actual fluids. 

## Running demos 

### Prerequisites

Ensure that you have python and a package manager (i.e. pip) installed. Requirements are listed in requirements, so once you're ready to install any packages, run the following in your terminal from the main directory: 

> pip install -r requirements.txt

### Scripts

There are two provided scripts. One shows a standing state of a few particles suspended solution. It can be run using the following command:

> python run_standing_display.py

The second simulation shows one particle being moved to another. It can be run using the following command:

> python run_meeting_display.py


## References

Stam, Jos. (2001). Stable Fluids. ACM SIGGRAPH 99. 1999. 10.1145/311535.311548. 

