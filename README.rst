=========
pykinetic
=========

------------------------------------------------------------------
A python library and command line tools for microkinetic modelling
------------------------------------------------------------------

.. project-description-start

This project consists of two command line apps "pykinetic-model.py" and 
"pykinetic-scan.py" to facilitate the generation of microkinetic models. Both 
apps generate a file that can be compiled/executed and that represents the 
kinetic model. It includes the "pykinetic" python library in case that either
the reactions or the chemical species are to be added programatically.

.. contents:: 
   :backlinks: none
   :depth: 2
   :local:

.. setup-instructions

Getting Started
---------------

These instructions will get you a copy of the project up and running on your
local machine.

Prerequisites
.............

- python >= 3.6 (3.7.3 or higher recommended)
- python library: setuptools
- python library: numpy (optional, needed to run the generated python scripts)
- python library: scipy (optional, needed to run the generated python scripts)
- c++ library: boost (optional, needed to run the generated c++ scripts)

   Note: proper documentation installation of the c++ library is being 
   currently being reviewed.

Installing the dependencies
...........................

python3.7 installation in Ubuntu 18.04 LTS

.. code:: shell-session

   $ sudo apt-get update
   $ sudo apt-get install python3.7 python3.7-dev


If for any reason it is not reachable:

.. code:: shell-session

   $ sudo add-apt-repository ppa:deadsnakes/ppa
   $ sudo apt-get update
   $ sudo apt-get install python3.7 python3.7-dev

Now you can skip the next step if you don't want to set up a virtual environment
the lines below correspond to python3.7 but can be also applied to python3.6 or 
higher (Remember to change "my_venv" for the actual path of the directory where 
you want the virtual environment)

.. code:: shell-session

   $ sudo apt-get install python3.7-venv
   $ python3.7 -m venv my_venv
   $ source my_venv/bin/activate

In the case of creating a virtual environment and running the 
previous "source" command the anytime you write "python" it is acting as alias 
for "python3.7" or whichever python executable you used for the virtual environment.
(To leave the virtual environment run the "deactivate" command at anytime)

If you have not created a virtual environment, all the following commands must 
explicitly use the python executable in which you are going to install the 
the package. 
Now we install the python default installer pip in our virtual environment

.. code:: shell-session

   $ python -m pip install pip
   $ python -m pip install --upgrade pip
   $ python -m pip install setuptools

If it proceeded without any errors (pip and setuptools should already be 
installed) otherwise please check the `pip`_ documentation

.. _pip: https://pip.pypa.io/en/stable/installing/

Now to install the optional dependencies (Skip if you do not want them):

.. code:: shell-session

   $ python -m pip install numpy scipy
   $ sudo apt-get install libboost-dev

Installing pykinetic
....................


Get the source code from github through git and install it.

.. code:: shell-session

   $ git clone https://github.com/maserasgroup-repo/pykinetic.git pykinetic
   $ python -m pip install pykinetic/

If you do not have git or do prefer to download manually the source 
code as a .zip or .tar.gz do it install it. 

.. code:: shell-session

   $ python -m pip install pykinetic-0.0.0.tar.gz

.. 
   
   Note: If you prefer to unpack it you can do it but it is not needed

Running the tests
.................

After installing you should be able to run the tests: 

.. code:: shell-session

   $ python -m unittest -v pykinetic.tests

Uninstalling pykinetic
......................

.. code:: shell-session

   $ python -m pip uninstall pykinetic


Developed with
--------------

- python 3.7.3
- Ubuntu 16.04 LTS and Ubuntu 18.04 LTS

.. examples-msg

Examples
--------

Please open the `Examples.rst <Examples.rst>`_ file in github to visualize the basic usage examples
or read the documentation.


.. project-author-license

Authors
-------

* **Raúl Pérez-Soto** - - https://github.com/rperezsoto


License
-------

(None currently)
