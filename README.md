----------
# Project
----------
**VENUS** development for etched-TJ and MJ-VCSELs.

Self-consistent 3D electro-opto-thermal solver for VCSEL simulations (a.k.a. VENUS). The original code dealt with oxide-confined VCSELs with some variations from the standard strcture, _i.e.,_ surface-relief VCSELs.

----------
How to use and brief code description
----------
Here there is a summary of how to launch a simulation in VENUS (either 1D and 3D). 

## Launch a simulation ##
- The full simulation must start with one of the `Main_*.m` scripts (please use some sort of significative name when creating new ones)
  - first the scripts asks the user to insert a string to append to the `.str` filename, needed to provide some significant name to the final result. All the results of a simulation are stored in the `out` folder; otherwise, the final result will be stored just with the structure filename (this can be useful in case of "not-significant" simulations, which therefore will be overwritten each time)
  - then, the script reads the set of parameters used for the simulation from a `settings_vari_*.m` file, containing among the others the name of strucutre, the type of simulation (electical, thermal, quantum, optical), the driving (both voltage and current driving (`mode.Idrive=1`) are allowed)
- The largest part of the code is inserted inside the script `SUB_drive`
  - the geometry of the structure is generated from the reading of the `*.str` file, and converted in a `geom` structure inside the script `Sub_str_from_VELM`, taht internally calls the scripts "StrRead*.m", each one for a different king of VCSEL
  - the geometry structure is then used as an input of the function `rectmesh.m`, where the multidimensional mesh is generated and stored in the `mesh` structure, together with a set of parameters
  -  the function `loadmacro.m` receives as input "mesh" and "mode" to infer to the nodes of the mesh the correct material properties 

### Maintainance ### 
- Pull this repo in a local repository, then create your own branch: name it as `[firt name letter]_[surname]_[idea of what the branch is needed for]`
  - Committing must be done by giving a proper name to the changes!
  - A **merge** request should be done only when some back-compatibility **tests** has been fulfilled, _i.e.,_ by running **at least** a VENUS quasi1D and 3D simulation on the structure with 1 TJ (this is simply done by running the `MAIN_C0.m`, therefore avoid touching that script and the related settings file `seetings_vari_C0.m`: please generate a copy and "play" with it!)
  - Please, make sure that the guy in charge of approving the merge request is aware of the request. Once the request is done, please avoid further modificying the branch scripts
- Please keep the code as similar as possible as its original version, _e.g.,_ the space before and after the `=` or spacing before the ending `;`; indent the `for` and `while` loops with `CTRL+E`

----------
Warnings and reminders
----------

- Please notice that in this version (20220620) the repo does NOT contain the folder "out/VELM", that is needed to the code for running and therefore must be created before launching any simulation
- The repository is set NOT to track the *.mat* files (important ones (<100 MB) could be added manually by using `git add -f [filename.mat]`): therefore, the LUT for the QW gain are not (and cannot) directly stored in the repo and should be taken from the shared folder [QW LUT](https://politoit-my.sharepoint.com/:f:/g/personal/s278303_studenti_polito_it/En5VLNrBVqxPrIt5yTuBJTwBZq80Eg-WhL5zbXJI_gsb-w?e=FC0zjL) (contact Alberto Gullino at alberto.gullino@polito.it), and placed into the folder "dati"
- Be careful: this code does NOT allow any space in the path, otherwise saving errors will be returned

----------
Acknowledgments
----------

The solver has been developped by Dr. Alberto Tibaldi and Dr. Pierluigi Debernardi from Politecnico di Torino (DET) and italian CNR, respectively. Since then, it has been extended by Dr. Marco Calciati and Ph.D. students Alberto Gullino, Valerio Torrelli and Martino D'Alessandro, all from DET department in Politecnico di Torino.
