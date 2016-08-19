Bacon age-depth modelling
--

This folder contain a script to introduce Bacon age-depth modelling. 
In this workflow, we pull geochronological data from Neotoma, create a Bacon input file, run Bacon, 
and then pass the created age-depth model back to our data object from neotoma.

To get Bacon to run, you need to have the script `Bacon.R` in your working directory. 
It must contain a folder names `Cores` and one called `bin`. 
If the version of `Bacon` in this folder does not work for you, you will need to get the 
version of Bacon you need from `http://www.chrono.qub.ac.uk/blaauw/bacon.html` and reorganize files as needed. 
This is annoying.
