# Mutualistic coevolution and community diversity favor persistence in metacommunities under environmental changes

This README_Cosmo_et_al_mutualistic_coevolution_metacommunities.txt file was generated on 2022-12-05 by Leandro Giacobelli Cosmo

GENERAL INFORMATION

1. Title of Dataset: 

Cosmo_et_al_mutualistic_coevolution_metacommunities

2. Author information:

Leandro G. Cosmo: Programa de Pós-Graduação em Ecologia, Departamento de Ecologia, Instituto de Biociências, Universidade de São Paulo - USP, São Paulo, SP, Brazil.

Lilian P. Sales: Departamento de Biologia Animal, Instituto de Biologia, Universidade Estadual de Campinas - UNICAMP, Campinas, SP, Brazil; Biology Department, Faculty of Arts and Science, Concordia University, Montreal, Canada.
 
Paulo R. Guimarães Jr.: Departamento de Ecologia, Instituto de Biociências, Universidade de São Paulo - USP, São Paulo, SP, Brazil.

Mathias M. Pires: Departamento de Biologia Animal, Instituto de Biologia, Universidade Estadual de Campinas - UNICAMP, Campinas, SP, Brazil.

Corresponding author: Leandro G. Cosmo, E-Mail: lgcosmo@usp.br

2. Information about funding sources that supported the collection of the data:

This study was financed in part by the Coordenação de Aperfeiçoamento de Pessoal de Nível Superior – Brasil (CAPES) – Finance Code 001. LGC thanks FAPESP for a PhD scholarship (São Paulo Research Foundation; grant # 2019/22146-3). LPS receives a Banting Postdoctoral Fellowship from the Government of Canada (Natural Sciences and Engineering Research Council; agreement # 01353). PRG is funded by CNPq (307134/2017-2), FAPESP (São Paulo Research Foundation; grant: # 2018/14809-0), and the Royal Society, London (CHL/R1/180156). MMP is funded by FAPESP (São Paulo Research Foundation; grant: # 2019/25478-7).

DATA & FILE OVERVIEW

1. File List: 

Data files:

main_simulation_results.csv

Scripts/Source functions:

main_functions.jl
aux_functions.jl
model_simulations.jl

DATA-SPECIFIC INFORMATION:

main_simulation_results.csv: full dataset containing the results of the numerical simulations used for the analyses in the main text. The variables in the dataset correspond to: 

(1) simulation - ID of the simulation
(2) species - ID of the species
(3) occupancy - Patch occupancy, defined as the proportion of the available patches occupied.
(4) local_tm - Average local trait matching with mutualistic partners
(5) local_em - Average local environmental matching.
(6) regional_tm - Average regional trait matching with all mutualistic populations.
(7) extinction_ratio - Average rate of extinction.
(8) colonization_ratio - Average rate of colonization.
(9) sp_extinct - Describes whether the species was extinct or not in the simulation. Cell value is equal to 1 if species went extinct and 0 otherwise.
(10) n_a - Number of animal species in the beginning of the simulation.
(11) n_p - Number of plant species in the beginning of the simulation.
(12) m - Strength of mutualisms as selective pressures.
(13) gene_flow - Fraction of gene flow among populations.
(14) alpha - Parameter that controls the shape of the trait matching function.
(15) gvar - Parameter that controls the additive genetic variance.
(16) rho - Parameter that controls the slope of species adaptive landscape.
(17) climchange - Amount of directional climate change in the simulation.

main_functions.jl/aux_functions.jl: Julia functions used to run the model numerical simulations.

model_simulations.jl: script to reproduce the numerical simulations of the model used in the main text.

USAGE INSTRUCTIONS:

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> Cosmo_et_al_mutualistic_coevolution_metacommunities

It is authored by Cosmo et al.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths. 

After installing everything, run the script "model_simulations.jl" located at the "scripts" folder.