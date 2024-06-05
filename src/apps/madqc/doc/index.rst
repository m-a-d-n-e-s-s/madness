.. maddft documentation master file, created by
   sphinx-quickstart on Tue Apr 30 14:32:26 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to maddft's documentation!
==================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation.rst
   quickstart.rst
   references.rst
   usage/calculation_structure.rst
   usage/qcengine_interface.rst
   usage/response_properties.rst

Multiresolution Quantum Chemistry Interface
-----------------------------------------------

The goal of this project is to provide a unified interface for computing with the various quantum chemistry applications available within MADNESS.
The motivation for this project is make MADNESS more accessible to the quantum chemistry community, providing them a platform to compute
benchmark quality results with minimal effort.  Towards easy access, we aim to make this interface compatible with QCEngine, given 
users the ability define and run their computations from python jupyter notebooks.

Currently, the interface is designed to work with :code:`moldft` and :code:`molresponse` in order to compute Hartree-Fock and DFT ground and response properties including excited states.
Future work will include additional quantum chemistry methods.

Why we need MADQCEngine?
------------------------

In the currently implemtation, all quantum chemistry methods are seperate executables, even though many of them depend on
each other.  For example, in order to compute a response property, which depends both on ground and response states, the current procedure
is to first compute the ground state using :code:`moldft` and then compute the response property using :code:`molresponse`.
To do this properly, a user necessarily needs to manage details such as
where to place the response calculation in relation the ground state calculation, and how to pass the ground state calculation.
Further complications arise when trying to compute response states at multiple frequencies, or higher order response properties which 
depend on a mixture of multiple ground and response states.  Currently the locations of all files where managed outside of the
code base which can be a cumbersome process.

Instead, MADQCEngine aims to provide a unified interface to manage all of these details for the user.  The user can simply define their
calculation, and the interface will manage the details of where to place the files in a consistent manner.  Additionally, the interface
will contain logic to manage the restart of files as well as whether or not to overwrite files during restarts.  For ease of use,
all outputs will be stored json files and will include a :code:`calc_paths.json` file which will contain the paths to all of the
files generated during the calculation.  


Getting Started
---------------

This documentation will guide you through the features, setup, and use of MADQCEngine. Here's how to get started:

1. **Installation**: Learn how to install MyProject on your machine. See :doc:`installation`.
2. **Quick Start**: Jump straight into using MyProject with our Quick Start guide. See :doc:`quickstart`.
3. **Reponse Properties** Learn how to compute HF/DFT response properties. See :doc:`usage/response_properties`.
4. **Calculation Structure** Learn how :code:`madqce` manages calculations and where to find the output files. See :doc:`usage/calculation_structure`.
5. **Interfacing with QCEngine** Learn how to interface with QCEngine. See :doc:`usage/qcengine_interface`.

Navigate through the documentation using the sidebar. Each section is designed to help you understand and use MyProject effectively.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`



