GUI Components
=================

The cBLUE GUI consists of a single window with a menu bar.

.. image:: ../images/GUI.png

Menu bar
********

The menu bar has the following categories.

.. csv-table:: cBLUE File Menu Options
	:header: Category, Options, Description
	:widths: 5, 6, 20
	
	File, Save Settings, Saves specified directories to a Json file
	.., Exit, Closes cBLUE		
	Sensor Model, Riegl VQ-880-G, "Currently, this is only a dummy menu option.  The senor model configuration for the Reigl VQ-880-G is hard-coded into cBLUE.  Development plans include refactoring the code to read sensor model information from a separate file and extending support to other lidar systems."
	Help, Documentation, Opens the cBLUE documentation in a web browser
	.., About, Displays the cBLUE splash screen

Data Directories
****************

cBLUE requires three directories to be set:

.. csv-table:: cBLUE File Menu Options
	:header: Directory, Description
	:widths: 6, 20
	
	Trajectory, contains the ASCII trajectory 
	LAS, contains the LAS files
	Output, where the output files will be created

.. _environ-label:
	
Environmental Parameters
************************

The subaqueous TPU calculations rely on characterizing two general environmental conditions:  water-surface wave conditions and turbidity.

Water Surface
-------------

The user has two options to characterize the water-surface. Both options use pre-computed coefficient lookup tables generated using Monte Carlo simulations.

* *Riegl VQ-880-G*
	This proof-of-concept option is currently limited in applicability. The lookup table used in this approach is generated from a water surface defined by a sampling of the first returns in a small area of the NOAA RSD Marcos Island project. Future development plans include deriving representative surface models for each flight line within each input tile.

* *Model (ECKV spectrum)*
	The lookup table used in this option is generated using a version of the Elfouhaily et al. (ECKV) directional gravity-capillary wave spectrum model that characterizes the distribution of wave frequencies using estimates of wind speed.
	
	.. csv-table:: ECKV Spectrum Options
		:header: options, wind speed range (kts)
		:widths: 10, 10
		
		Calm-light air, 0-2
		Light Breeze, 3-6
		Gentle Breeze, 7-10
		Moderate Breeze, 11-15
		Fresh Breeze, 16-20

Turbidity
---------

The user can chose among 5 turbidity classes.

.. csv-table:: Turbidity (Kd_490) Classes
	:header: class, value range (:math:`m^{-1}`)
	:widths: 10, 10
	
	Clear, 0.06-0.10
	Clear-Moderate, 0.11-0.17
	Moderate, 0.18-0.25
	Moderate-High, 0.26-0.32
	High, 0.33-0.36

VDatum Region
*************

The user has the option to select a VDatum region, to include the corresponding maximum cumulative uncertainty (MCU_) in the final TPU calculations.  

.. _MCU: https://vdatum.noaa.gov/docs/est_uncertainties.html

Process Buttons
***************

	.. csv-table:: Process Buttons
		:header: Button, Function
		:widths: 20, 50
		
		Load Trajectory File(s), Loads the trajectory data from the ASCII sbet files located in the user-specified Trajectory directory into a single Pandas dataframe
		Process TPU, Initiates the TPU calculations and creates Las files with calculated TPU as variable length record (VLR) extra bytes

