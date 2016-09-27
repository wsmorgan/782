1D Quantum Pontentials
======================

The basis solution requires knowledge of the potential, defined in a
potential config file.

Configure Sections
------------------

A potential is defined using the section `[parameters]` and
`[regions]`:

- **[parameters]** has a list of *lowercas* parameter names as
  options with default values.
- **[regions]** options are region numbers, values have the form:

  .. code-block:: bash

    start, stop | region value (or function)

  `start` and `stop` can be any valid python code where the variable
  names should be defined in `[parameters]`. The value or function can
  also refer to any variables; for functions use th `lambda` syntax:

  .. code-block:: python

     lambda x: v0*(x-shift)**2 - numpy.exp(x)

  Functions can use any `numpy` or `operator` modules in addition to the
  parameters.

API Documentation
-----------------

.. automodule:: basis.potential
   :synopsis: How to configure a 1D potential using the configuration file.
   :members:
