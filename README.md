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
<br>
- File: example_fiberLengthsCal.m <br>
  Objective: Calibrates muscle fiber lengths with provided digitalized values.
  Outcome: Optimized optimal fiber lengths, tendon slack lengths, and tendon stiffness
 
- File: example_passMomCal.m <br>
  Objective: Calibrates passive angle-moment relationships at various ankle, knee, and hip joint angles. Digitalized data provided.
  Outcome: Optimized muscle passive force-length curves
  
- File: example_fiberLength_PassCal.m <br>
  Objective: To calibrate muscle fiber lengths and passive angle-moment relationships. Digitalized data provided.
  Outcome: Optimized optimal fiber lengths, tendon slack lengths, tendon stiffness, and muscle passive force-length curves
  
- File: example_exoskeleton.m <br>
  Objective: To compute the optimal assistance based on spring-like and motor-like actuation systems to support plantarflexors. Additionally, we added an example where you can create your own assistive moment to support knee extensors. You need to run "example_fiberLength_PassCal.m" as this example requires optimized values previously computed.
  Outcome: Optimal assistive moment and corresponding muscle-tendon dynamics.
  
<br>
Also, each script provides a description of the computation performed.
<br>

#### About the experimental data 
We provide a scaled musculoskeletal model and motion data: inverse kinematics, inverse dynamics, of a subject walking at a preferred speed. Digitalized data is also provided. Please take a look at our studies for more details.
<br>

## Contact 
For further information or questions, feel free to reach me by email: ailp@kth.se
