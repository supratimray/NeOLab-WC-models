# NeOLab-WC-models
Codes used to organize, specify and simulate Wilson-Cowan type models with different inputs
This repository is made of stripped-down versions of codes used to run WC models in Krishnakumaran et al 2022 and Krishnakumaran and Ray 2023.

## Get Started
Add all folders to path and Run runexamples.m

----------------------------------------------------------------
# NeoLab-WC-models-basic is another repository which implements some simple WC models and displays the responses for constant inputs given to the excitatory and inhibitory populations. This repository uses a more advanced coding structure which does the following:
1. Different models are incorporated in a single structure. This is done by first creating generic classes(see Neurons.m and GenerateNeuronPopulation.m under MAP/DataStructure and backend MAP/Simulation - Network Classes). Subsequently, different "model cards" under MAP/ISN model library MAP are configured.

2. Allows the inputs to be noisy.

3. Allows time varying inputs to E and I populations.
4. Can potentially be expanded to other models as well, such as SSN.

# Contents

## DataStructure and backend MAP
- Defines the model object - data structure and parametrization codes
- Defines simulation backend codes adaptable to different WC models tested to allow defining them under a common framework

## ISN model library MAP
- Contains "WC Model Cards" - functions with previously used WC model definitions that return configured model object - summarized in parameter tables in papers
- Output objects need only be configured a bit further to specify simulation time and few more details.

## Example codes - ISN model scripts
- Some example simulations used previously for different WC model and functions used therein

## misc-PublishableFigureFormatting
- contains code used in some example codes to auto-format figure axes and plots/scatters. 

----------------------------------------------------------------

### Some notes on Output structure

Outputs will be stored in newer folder (with model name unless overridden in run codes)

The simulation model object containing all the parameters values for each simulation iteration, definitions of non-linear activation functions or input update functions and input file information, if any, along with the simulated output timeseries can be found in the output of these example codes as popln or JS_pop variables

Before loading these files, it is required to add "DataStructure and backend MAP" directory to path so that Matlab recognizes the objects and parses them meaningfully.

The object could also be used for further extended simulations, in place, by redefining its tspan and relevant input information.




