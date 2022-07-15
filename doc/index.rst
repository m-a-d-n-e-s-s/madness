.. madness documentation master file, created by
   sphinx-quickstart on Tue Jul 12 21:54:04 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to MADNESS
==================
Multiresolution Adaptive Numerical Environment for Scientific Simulation

MADNESS provides a high-level environment for the solution of integral and differential equations in many dimensions using adaptive, fast methods with guaranteed precision based on multi-resolution analysis and novel separated representations. There are three main components to MADNESS. At the lowest level is a new petascale parallel programming environment that increases programmer productivity and code performance/scalability while maintaining backward compatibility with current programming tools such as MPI and Global Arrays. The numerical capabilities built upon the parallel tools provide a high-level environment for composing and solving numerical problems in many (1-6+) dimensions. Finally, built upon the numerical tools are new applications with initial focus upon chemistry, atomic and molecular physics, material science, and nuclear structure.



.. toctree::
   :maxdepth: 2

   runtime
   numerical_library
   quantum
   use MADNESS

..
  - :doc:`runtime`
  - :doc:`numerical_library`
  - :doc:`quantum`
  - :doc:`use_madness`

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
