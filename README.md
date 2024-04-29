# The State of Goa energy system model
![Alt text](./figures/F1_methodology.png)
*Figure 1. State of Goa-model development methodology*

This repository presents a open-source energy models explicitly developed for the State of Goa (India) Kenya. The energy model consists of several sectors such as, building (residential + commercial), cooking, industry, agriculture, fisheries, and transport. This energy system model was develped to be used wiht OSeMOSYS-pulp [^1].

## Input data
The input data is a **.xlsx* file located in:
```
./model/Input_Data/
```
The input data contains all the parameters necessary to run the optimization. Documentation regarding the parameters are similar to those found in the original [OSeMOSYS-GNU](https://osemosys.readthedocs.io/en/latest/) [^2].

## Model file
The model file is a updated version of OSeMOSYS pulp. Several updted were made to speed up optimization regarding matrix gerenation time, and post-processing. Now, this version is comparable to the short version of OSeMOSYS-GNU. The model file is located in:
```
./model/OSeMOSYS.py/
```
Other dependencies and functions are located in:
```
./model/utils/
```

## Run the model
The working directory should be set in ./model. It is advisible to create a new enviroment and install the libraries such as pulp.
```
conda create --name pulp python=3.8
```
Activate the new environment to use it
```
conda activate pulp
```
 run the following code, and to use other solvers type instead of cplex, gurobi, or cbc. The result file is a **.csv* file save in the path *./Output_Data/GOA_COMPLETE_updated_results.csv*
```
python OSeMOSYS.py -i GOA_COMPLETE_updated.xlsx -s cplex -o csv
```

## Vizualize results
In this version of the model, a python notebook is provided to vizualize results. The notebook is located in the following path:
```
./scripts/visualizarion_csv_Goa.ipynb
```
[^1]: D. Dreier, M. Howells, Osemosys-Pulp: A Stochastic Modeling Framework for long-term Energy Systems Modeling, Energies. 12 (2019) 1382. doi:10.3390/en12071382. 
[^2]: M. Howells, H. Rogner, N. Strachan, C. Heaps, H. Huntington, S. Kypreos, et al., OSeMOSYS: The Open Source Energy Modeling System, Energy Policy. 39 (2011) 5850–5870. doi:10.1016/j.enpol.2011.06.033. 