# The State of Goa energy system model
![Alt text](./figures/F1_methodology.png)
*State of Goa-model development methodology*

This repository presents a open-source energy models explicitly developed for the State of Goa (India) Kenya. The energy model consists of several sectors such as, building (residential + commercial), cooking, industry, agriculture, fisheries, and transport. 

 The first model is the [Whole Energy System for Kenya (WESM)](https://github.com/ClimateCompatibleGrowth/osemosys_kenya), encompassing power generation, industrial, and transportation sectors, among others. The second model is a Climate, Land, Energy, and Water system [(CLEWs)](https://github.com/robertodawid/Kenya_Clews/tree/main) model, focusing on intricate interdependencies between the energy and land systems within the Kenyan context. The CLEWs model presents the interlinkages evident in sectors such as cooking (representing energy) and agriculture (representing land). While both models address aspects of the energy system, the CLEWs model uniquely integrates additional dimensions such as land and water systems.
Merging these models has the potential to capture additional interactions between various systems. For instance, the utilization of fossil fuels extends beyond energy generation to activities in the land and water domains, including mechanization and pumping, respectively. Similarly, electricity finds application in residential settings for irrigation and public water distribution, showcasing interdependencies across sectors.

## Read and Print the results in easy names.

This short Jupyter notebook contains the script to read the results after using a least-cost optimization approach. The Pulp version of OsEMOSYS [^fn]  was used to read, create the lp file and solve the energy model for the State of Goa.

[^fn]: D. Dreier and M. Howells, “OSeMOSYS-PuLP: A Stochastic Modeling Framework for Long-Term Energy Systems Modeling,” Energies, vol. 12, no. 7, p. 1382, Apr. 2019, doi: 10.3390/en12071382.
