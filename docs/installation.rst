Installation
============

Setting up the environment
--------------------------

First, clone this repo locally:

.. code-block:: bash
   
   git clone https://github.com/spestana/goes-ortho
   cd goes-ortho


Using `conda <https://docs.conda.io/projects/conda/en/latest/index.html>`_ or `mamba <https://mamba.readthedocs.io/en/latest/>`_

.. code-block:: bash

   conda env create -f environment.yml
   conda activate goesenv
   pip install -e .

If you are using Jupyter Notebooks, you may also need to run the following to use this environment within your notebooks:

.. code-block:: bash

   ipython kernel install --user --name goesenv

----
