# Exoskeletons & Experiment-guided parameter calibration
## Overview
In this repository, we provide a simulation framework to calibrate muscle-tendon parameters based on prior recorded data and simulate muscle-tendon dynamics based on the assistive moment, i.e., ideal exoskeleton assistance. This  framework uses trajectory optimization to simulate muscle dynamics and solve the muscle redundancy. Further information about our study and computational methods can be found in the publications:
- Experiment-guided calibration of muscle fiber lengths and passive forces
- Springs vs. motors: Optimal assistance to support walking across speeds
<br>
This is a work in collaboration with VU Amsterdam. The MatLab scripts are derived from the muscle redundancy solver from KULeuvenNeuromechanics.
<br>

## Installation instructions
You need to use Matlab to run the scripts and download three software packages:
* Install Opensim. Our algorithms have been tested with OpenSim 4.1, an open-source software for musculoskeletal modeling and simulations. See the following link for installation: https://simtk.org/frs/?group_id=91 
* Add the tool “Scripting with Matlab.” Get access to the functions of OpenSim from Matlab. See the following link for installation: https://simtk-confluence.stanford.edu:8443/display/OpenSim/Scripting+with+Matlab
* Install Casadi. Open-source tool for nonlinear optimization and algorithm differentiation. See the following link for installation: https://web.casadi.org/ 
<br>
This code is based on the Muscle Redundancy Solver: https://github.com/KULeuvenNeuromechanics/MuscleRedundancySolver. It might be a good idea first to run the example they provided to verify that your installation is correct, yet it is optional.
<br>

## Simulation framework
Our simulation framework calibrates muscle-tendon parameters and simulates optimal assistive moment and its corresponding muscle-tendon dynamics interaction. Each tool can be used independently. The inputs are a scaled musculoskeletal model computed from OpenSim and additional inverse kinematics and dynamics solutions, depending on the simulation objectives.
<br>

### Experiment-guided muscle-tendon parameter calibration
This tool optimizes muscle-tendon parameters to better resemble in vivo observations of muscle fiber lengths and angle-moment relationships reported in the literature. It calibrates optimal fiber length, tendon slack length, and tendon stiffness to match reported digitalized images from ultrasound and muscle passive curves to match reported in vivo experimental angle-moment relationships. 

### Simulation of the optimal assistive moment and muscle-tendon interaction
This tool optimizes an assistive moment that minimizes muscle activations. You can specify the assistive moment based on spring-like or motor-like actuation at various muscle groups within the musculoskeletal model.
<br>

 ## Run example
The scripts can be adapted to your research. We facilitate an example and experimental data for each script; see the folder “Code and example.” To run any of the tools, you first must follow the “Installation instructions.” Also, verify that Casadi and the main folder (Code and example) paths are added to your current Matlab session.
- Example muscle-tendon parameter calibration
Run "muscle-"


The three scripts combined allow you to compute the lower limb’s metabolic rate estimations using a simulation workflow and metabolic energy model used in our study. In order to obtain such a result, you must follow the following steps:
2) Run the script “Calibration_passiveForces.” This script serves to calibrate passive force parameters in a generic musculoskeletal model. As a result, it will provide a summary of the calibrated parameters, store the results, and plot the graph: experimental, calibrated, and generic passive torque-angle curves.
3) Run the script “Pipeline_simulationFramework.” This script serves to select the simulation workflow you want to use. It specifies the features for the simulation and uses the function “setupAndRun” to run the code. As a result, it will provide multiple graphs with information about the muscle activations, normalized fiber length, fiber length, and reserve actuators. Note: You must run “Calibration_passiveForces” first, as this script requires generic or calibrated passive parameters.
4) Run the script “Compute_metabolicRates.” This script serves to compute metabolic rates based on the muscle excitation, states, and state derivatives previously obtained from the MRS. We provide the implementation of various metabolic energy models; you can select any. As a result, it will provide estimates of the metabolic rates, work rates, and heat rates per muscle. Note: You need to run “Pipeline_simulationFramework,” as this script requires muscle-tendon states to compute the metabolic rates.
<br>
Also, each script provides a description of the computation performed.
<br>

#### About the experimental data 
We provide a scaled musculoskeletal model and motion data: inverse kinematics, inverse dynamics, and EMGs of a subject walking at slow, normal, and fast speeds.  Please take a look at our studies for more details.
<br>

## Contact 
For further information or questions, feel free to reach me by email: ailp@kth.se
