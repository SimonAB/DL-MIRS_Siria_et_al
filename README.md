[![DOI](https://zenodo.org/badge/296672842.svg)](https://zenodo.org/badge/latestdoi/296672842)

# Rapid age-grading and species identification of natural mosquitoes for malaria surveillance

Doreen J. Siria, Roger Sanou, Joshua Mitton, Emmanuel P. Mwanga, Abdoulaye Niang, Issiaka Sare, Paul C.D. Johnson, Geraldine Foster, Adrien M.G. Belem, Klaas Wynne, Roderick Murray-Smith, Heather M. Ferguson, Mario González-Jiménez, Simon A. Babayan, Abdoulaye Diabaté, Fredros O. Okumu, and Francesco Baldini

**The malaria parasite, which is transmitted by several
  _Anopheles_ mosquito species, requires more time to reach its
  human-transmissible stage than the average lifespan of mosquito vectors.
  Monitoring the species-specific age structure of mosquito populations is
  critical to evaluating the impact of vector control interventions on malaria
  risk. We developed a rapid, cost-effective surveillance method based on deep
  learning of mid-infrared spectra of mosquito cuticle that simultaneously
  identifies the species and age class of three main malaria vectors in natural
  populations. Using spectra from over 40,000 ecologically and genetically
  diverse _An. gambiae_, _An. arabiensis_, and _An. coluzzii_
  females, we developed a deep transfer learning model that learned and
  predicted the age of new wild populations in Tanzania and Burkina Faso with
  minimal sampling effort. Additionally, the model was able to detect the
  impact of simulated control interventions on mosquito populations, measured as
  a shift in their age structures. In the future, we anticipate our method can
  be applied to other arthropod vector-borne diseases.**

## Resources

- [Code for converting spectra](code/Spectra%20conversion%20script/Bad%20Blood%201.2.ipynb)
- [Code for UMAP Clustering](code/UMAP/UMAP_DL-MIRS_data.ipynb)
- [Code for convolutional neural net](code/CNN/)
- [Code for power analyses](code/Power%20analyses/MIRS_AgeStructure_PowerSim_v05.R)
- [Code and instructions for calculating age structure proportions based on gonotrophic cycles](code/AgeStructureProportions)

## General system requirements

Developed on:

- Operating systems: macOS 10-12; Windows 7; Linux Ubuntu
- Hardware: CPU Intel Core i7, 16-64 GB RAM
- Specialised hardware for deep learning: GPU - TITAN Xp 12GB

## Installation guide

- Requires Python 3.x (we recommend [Anaconda](https://www.anaconda.com/products/individual)) and [R v. 4 and above](https://cran.r-project.org), although older versions are likely to work too. Details of packages used are provided in the accompanying paper's methods section and in README files where relevant within each of the folders listed under [resources](#Resources)
- Download this repository to your computer
- Download the data from Enlighten (link coming soon)
- Place the data in the same directory as the scripts

## Instructions for use

Scripts for processing MIRS data, training models, and reproducing are provided as jupyter notebooks. Before you can run the scripts as provided, take care to start by updating the file paths pointing to the source data as per your directory structure.


Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
