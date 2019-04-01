.. _subaqu-label:

Subaqueous Component
********************

The subaqueous component stochastically computes the TPU of the laser pulse on the seafloor using Monte Carlo ray tracing based on surface modeling and estimates of scattering and absorption (see the figure below).  To minimize the computational complexity of performing Monte Carlo simmulations for large numbers of data points (> 1 billion), the algorithm relies on pre-computed coefficient lookup tables.  

The figure below illustrates the general concept of the subaqueous TPU calculations.  Given assumed sea state (and turbidity) conditions, large numbers of light rays are propagated, or traced, through the water column.  The spread of the resulting distribution of depths provides an estimate of uncertainty.

Please refer to the `cBLUE Monte Carlo Manual`_ for detailed information about the subaqueous Monte Carlo simulations and corresponding lookup tables.

.. _`cBLUE Monte Carlo Manual`: ../html/_static/MonteCarlomanual.pdf

.. image:: ../images/SubaqueousMonteCarlo.gif

Image Credit: Firat Eren
