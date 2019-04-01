.. image:: ./cBLUE_readme.gif

comprehensive Bathymetric Lidar Uncertainty Estimator (cBLUE)
=============================================================

cBLUE is a tool to calculate the total propagated uncertainty of bathymetric lidar data.  Currently, the tool only supports the Riegl VQ-880-G system, but development plans include adding support for other lidar systems.  The theoretical foundation for cBLUE is documented in Eren et al. (in press) [#]_.

cBLUE is a collaboration among the `Remote Sensing Division`_ within `NOAA's National Geodetic Survey`_ and researchers within the School of Construction and Civil Engineering at Oregon State University (`Parrish Research Group`_) and the Center for Coastal and Ocean Mapping/Joint Hydrographic Center (`CCOM/JHC`_) at the University of New Hampshire. 

.. _`NOAA's National Geodetic Survey`:  https://www.ngs.noaa.gov

.. _`Remote Sensing Division`:  https://www.ngs.noaa.gov/RSD/rsd_home.shtml

.. _`Parrish Research Group`: http://research.engr.oregonstate.edu/parrish/

.. _`CCOM/JHC`: http://ccom.unh.edu/about-ccomjhc

Installation
============

cBLUE, created with Python 3, is currently provided as a collection of Python modules (.py files) and supporting text files.  The cBLUE GUI is accessed by running the CBlueApp.py file via the command line or a Python IDE. 

Dependencies
------------

The recommended way to ensure that all of the necessary dependencies are loaded is to create a conda environment from the text file cBLUE_install.txt using the following command, 

	conda env create --file <evn file name>
	
where <env file name> is the path to cBLUE_install.txt, which is included in the cBLUE GitHub release.

Starting the GUI
----------------

The GUI can be initiated from the command line with the following command (specify the full path of CBlueApp.py if the current directory is not the location of CBlueApp.py):

NOTE: If there are multiple Python installations on the computer running the GUI, ensure the desired installation is specified in the path environment variables or explicitly specify the full path of the relevant python executable in the command. For example:

	(command line)> C:\\Python32\\64-bit\\python.exe CBlueApp.py

The GUI can also be initiated by running the CBlueApp.py file from a Python IDE, such as IDLE or PyCharm.

Documentation
=============

Access cBLUE documentation (under construction) by opening docs/index.html (from your local repo) in a web browser.  

.. rubric:: Footnotes

.. [#] Eren, F., Jung, J., Parrish, C. E., Forfinski-Sarkozi, N., and Calder, B. R., “Total Vertical Uncertainty (TVU) modeling for topo-bathymetric lidar systems”, American Society for Photogrammetry and Remote Sensing (ASPRS). In Press.
