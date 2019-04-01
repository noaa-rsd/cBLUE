Algorithm Workflow
------------------

The cBLUE algorithm workflow consists of three main steps:

1. Form Sensor Model Observation Equation

	A laser geolocation equation is formed based on characteristics of the specified lidar sensor model.  The symbolic math library *SymPy* (https://www.sympy.org/en/index.html) is used to construct the laser geolocation equation from specified laser and airplane rotation matrices.

2. Generate Jacobian

	The general equation of the Jacobian, or the matrix of partial derivatives, of the laser geolocation equation is calculated using SymPy.

3. Propagate Uncertainty

	Once the general equation of the Jacobian is calculated, the uncertainties of the component variables are propagated for each data point, per flight line, per Las tile, using the following steps:

	* Merge the Las data and trajectory data (Merge class)
	* Calculate :ref:`subaer-label` (Subaerial class)
	* Calculate :ref:`subaqu-label` (Subaqueous class)
	* Combine subaerial and subaqueous TPU (Tpu class)
	* Export TPU as Las extra bytes (Tpu Class)
	