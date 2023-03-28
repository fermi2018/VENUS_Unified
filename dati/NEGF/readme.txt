In this folder are contained the interpolated results coming from NEGF simulations used to investigate the tunneling current across the TJ. 

To generate them, move to folder: 20211120_NEGFsimulations. There, a set of script must be run.

- Main_1_str_NEGF_TJ.m: 
	- read the str file from which a cut of the TJ is obtained ("strName"). In principle, this is the same file used in VENUS (copy and paste it in subfolder D1ANA/generageom). However, some adjustments are required: refer to "btj_GaAs.str" as a working file!
	- set the bias range to investigate with NEGF
		- reverse bias (where tunneling is relevant) is applied with a negative sign. The saved workspace is used for the interpolation; other files are saved to visualize NEGF quantities as outputs
		- forward bias is applied with a positive sign (typical pn junction operation)
	- while the simulation runs, the band diagram should be visualized: be careful that the bands overlook. At the end of the simulation, J_NDD (DD current) should be very close to zero!
	
- Main_2_Visualize_NEGF_Solutions.m:
	- set again the same str name as input ("strName")
	- shows GBTBT from NEGF. Other quantities (LDOS, spectral curr. dens., ...) can be visualized with inner functions
	
- Main_3_Plot_Tunneling_IV_and_Generate_Interp_From_State.m
	- set again the same str name as input ("strName")
	- load the Workspace from Main_1_str_NEGF_TJ to recover the bias window
	- at each bias point loads a ['State_TJ_',strName,'Vbias=.mat'] file containing the NEGF results at each bias
	- extracts: I_tunneling (current density), GBTBT (tunneling rate), VTJ (voltage drop across the TJ), node (TJ nodes)
	- save them in an ['Interp_',strName,'.mat'] file
	- plots JV and IV (using mode.Area) from NEGF
	
- Main_3_1_Plot_Tunneling_IV_and_IVfit_LOG_AG.m
	- set again the same str name as input ("strName")
	- pDeg: fitting order of the logarithmic fit
	- clog: contains the fitting coefficients
	- I0: current to be subtracted to the fitting at 0 V (current must be zero at 0 V!)
	- save them in an ['Interp_',strName,'_log.mat'] file, USED IN THE FULL DD SOLVER (VENUS/D1ANA)
	
Eventually, in VENUS, ['Interp_',strName,'_log.mat'] loaded in LUT_GBTBT.m to store the fitting coefficients in mode structure
	- use TJcharacteristic_test.m to see the IV from NEGF!