Installation
============

The latest release of cBLUE can be found at https://github.com/forkozi/cBLUE/releases.

cBLUE, designed to run on 64-bit Python 3, is currently provided as a collection of Python modules (.py files) and supporting text files.  The cBLUE GUI is accessed by running the CBlueApp.py file via the command line or a Python IDE. 

cBLUE Files
-----------

The following table lists the main classes comprising cBLUE:

.. csv-table:: Main cBLUE classes
	:header: file, description
	:widths: 14, 30

	CBlueApp.py, defines and initiates the GUI
	Sbet.py, loads the ASCII trajectory files (or "sbets")
	Las.py, loads the las files
	Merge.py, merges the trajectory and las data based on timestamps
	SensorModel.py, defines and gives access to lidar sensor model
	Jacobian.py, forms and evaluates the Jacobian of a sensor model's laser geolocation equation
	Tpu.py, coordinates the TPU workflow
	Subaerial.py, calculates the subaerial TPU
	Subaqueous.py, calculates the subaqueous TPU
	Datum.py, loads the zone-specific VDatum zone MCU

Dependencies
------------
The recommended way to ensure that all of the necessary dependencies are loaded is to create a conda environment from the text file cBLUE_install.txt using the following command, 

	conda env create --file <evn file name>
	
where <env file name> is the path to cBLUE_install.txt, which is included in the cBLUE GitHub repository.

.. warning::

	Due to a current issue with the most recent version of conda (https://github.com/conda/conda/issues/8404), the above method will not install the pip packages listed in the .yml file.  To install the necessary pip packages, use :code:`pip install -r cBLUE_env_pip.txt` in addition to the :code:`conda env create --file cBLUE_env.yml` command listed above.

Major dependencies are summarized in the table below:

.. csv-table:: Main Dependencies
	:header: package, comment
	:widths: 10, 30

	numpy, most of the matrix calculations are based on numpy arrays
	numexpr, used to accelerate numpy-array calculations of large arrays
	sympy, used to symbolically form geolocation equation and form/evaluate Jacobian
	pandas, certain operations are based on Pandas dataframes
	pathos, multiprocess framework (if user specifies multiprocessing)
	laspy, used to read LAS files
	Tkinter, used to create GUI

Starting the GUI
----------------
The GUI can be initiated from the command line with the following command (specify the full path of CBlueApp.py if the current directory is not the location of CBlueApp.py):

NOTE: If there are multiple Python installations on the computer running the GUI, ensure the desired installation is specified in the path environment variables or explicitly specify the full path of the relevant python executable in the command. For example:

	(command line)> C:\\Python27\\64-bit\\python.exe CBlueApp.py

The GUI can also be initiated by running the CBlueApp.py file from a Python IDE, such as IDLE, PyCharm, or Microsoft Visual Studio.

Development plans include packaging all of the necessary files into a single-file executable.
