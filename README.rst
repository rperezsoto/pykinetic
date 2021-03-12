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
(Remember to change "my_venv" for the actual path of the directory where you
want the virtual environment)

.. code:: shell-session

   $ sudo apt-get install python3.7-venv
   $ python3.7 -m venv my_venv
   $ source my_venv/bin/activate

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


Get the source code either git or download and unpack it into the "pykinetic"
folder.

.. code:: shell-session

   $ git clone https://github.com/maserasgroup-repo/pykinetic.git pykinetic

Now install pykinetic

.. code:: shell-session

   $ python -m pip install pykinetic/


Installing with the -e option before pykinetic will make that all the changes in
the source files will have have effect when you call them through their alias.
However, you have to manually clean the folder generated in case of uninstalling
the package.

Running the tests
.................

After installing the simplest way to run the tests is go to the tests folder and
run: 

.. code:: shell-session

   $ python -m unittest -v test_*.py

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
