=========
pykinetic
=========

------------------------------------------------------------------
A python library and command line tools for microkinetic modelling
------------------------------------------------------------------

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.8053050.svg
   :target: https://doi.org/10.5281/zenodo.8053050

.. project-description-start

Pykinetic is a project that aims to facilitate the automation and generation 
of kinetic and microkinetic models. Pykinetic acts as an intermediate layer 
between the user and the mathematical model that is then simulated to obtain the
kinetic data. This functionality can be done either through the use of the 
pykinetic python library or through the use of the command line apps: 
"pykinetic-model.py" and "pykinetic-scan.py". The minimal source code to run the 
simulation is the output of pykinetic and currently python and c++ are supported. 
Compiling and executing that code leads to the simulation of the chemical system. 

.. project-description-end

.. contents:: 
   :backlinks: none
   :depth: 2
   :local:


Getting Started
---------------

These instructions will get you a copy of the project up and running on your
local machine.

.. setup-instructions-start

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
previous "source" command, anytime you write "python" it will act as an alias 
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

Pykinetic can be directly installed through pip: 

.. code:: shell-session

   $ python -m pip install pykinetic

Get the source code from github through git and install it.

.. code:: shell-session

   $ git clone https://github.com/maserasgroup-repo/pykinetic.git pykinetic
   $ python -m pip install pykinetic/

If you do not have git or do prefer to download manually the source 
code as a .zip or .tar.gz do it install it. 

.. code:: shell-session

   $ python -m pip install pykinetic-0.1.0.tar.gz

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

.. setup-instructions-end

Developed with
--------------

- python 3.7.3
- Ubuntu 16.04 LTS and Ubuntu 18.04 LTS


Examples
--------

Please open the `Examples.rst <Examples.rst>`_ file in github to visualize the 
basic usage examples or read the documentation.

Authors
-------

.. project-authors-start

List of main developers and contact emails:  

*  Raúl Pérez-Soto [
   `ORCID <https://orcid.org/0000-0002-6237-2155>`__ ,
   `Github <https://github.com/rperezsoto>`__ ,
   `email <rperezsoto.research@gmail.com>`__ ]
*  Sergio Pablo-García [
   `ORCID <https://orcid.org/0000-0002-3327-9285>`__ , 
   `Github <https://github.com/spgarcica>`__ , 
   `email <spgarcica@gmail.com>`__]
*  María Besora [
   `ORCID <http://orcid.org/0000-0002-6656-5827>`__ ,
   `Github <https://github.com/BesoraMaria>`__ ,
   `email <maria.besora@urv.cat>`__ ] 
*  Feliu Maseras [
   `ORCID <http://orcid.org/0000-0001-8806-2019>`__ ,
   `Github <https://github.com/maserasgroup-repo>`__ , 
   `email <fmaseras@iciq.es>`__] 

.. project-authors-end

License
-------

.. project-license-start

pykinetic is freely available under an `MIT <https://opensource.org/licenses/MIT>`__ License

.. project-license-end
