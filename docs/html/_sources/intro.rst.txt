Overview
========

.. toctree::
	:maxdepth: 1 

This tool computes the vertical total propagated uncertainty (TPU) of bathymetry acquired with a Riegl VQ-880-G topobathymetric lidar (other lidar systems will be included in future versions).  The algorithm consists of subaerial and subaqueous components (see the figure below).  Each component lends itself to a different approach to uncertainty propagation because of the relative complexity of the factors influencing the laser pulse travel path.  Whereas the subaerial portion is a well-defined geometric problem that can be addressed using standard geomatics techniques, the subaqueous portion uses a Monte Carlo approach  to model the complex interactions of light with water that are difficult to model analytically. 

.. toctree::

	algorithm_overview
	subaer
	subaque

.. image:: ../images/GeneralApproach.png

Image Credit: Chris Parrish
