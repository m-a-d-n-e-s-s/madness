=======================
Properties
=======================


Properties can be computed with the by selecting driver as property.

.. code-block:: yaml

  task
  driver "property"
  method "tddft"

In the `property` block one defines the properties to be computed as well as methods controlling the computation.
For example, to compute the frequency-dependent polarizability one would write:


.. code-block:: yaml

  property
    polarizability
      tddft # controls the tddft linear response solver
      frequency_range 0.0 0.1 10


