.. coloring-start 

.. raw:: html

   <style> .red {color:#c76160; font-weight:bold; font-size:16px} </style>
   <style> .green {color:#3d642d; font-weight:bold; font-size:16px} </style>

.. role:: red

.. role:: green

.. coloring-end

==================
pykinetic examples
==================

.. contents:: Table of Contents
   :backlinks: none
   :local:


File Formats
------------

Compounds
.........

.. compounds-start

General format:
***************

Each line of the compounds file follows the format: ::

   # Line that starts with "#" is considered as comment and ignored
   compound_label    energy    unit    scan

* The label and the energy columns are always required. 
* The unit may be specified as command line argument. 
* The "scan" keyword allows to mark a certain compound to produce several models while changing the energy of everything that is marked by the scan keyword. 
* The number of empty spaces or tabs is not strict as long as all the columns are separated by at least one (tab or empty space).

Compound label rules:
*********************

The compound's label is relatively free. The only exceptions are those that 
include any symbol used in the reactions file to specify the reactions:

* must not have emtpy spaces within the label.
* "!" symbol must not be used.
* "+" may be included as long as it is not surrounded by empty spaces or tabs.
* Any reaction symbol surrounded by empty spaces or tabs.

+---------------------------+---------------------------+
|   **Accepted Labels**     |   **Non-Valid Labels**    |
+===========================+===========================+
| :green:`[H2O···H2O]`      | :red:`[H2O H2O]`          |
+---------------------------+---------------------------+
| :green:`H2O+H2O`          | :red:`H2O + H2O`          |
+---------------------------+---------------------------+
| :green:`[LCu(OTf)2]Na+`   | :red:`[LCu(OTf)2]Na +`    |
+---------------------------+---------------------------+
| :green:`Enz=>inactive`    | :red:`Enz => inactive`    |
+---------------------------+---------------------------+
| :green:`花火`             |          :red:`!`         |
+---------------------------+---------------------------+

Compounds file example:
***********************

.. code:: none

   A                0.0 kcal/mol
   B                0.0 kcal/mol
   Cat              0.0 kcal/mol
   [A···B]          2.0 kcal/mol
   [A···B···Cat]    4.0 kcal/mol
   [C···Cat]        3.0 kcal/mol
   C               -2.0 kcal/mol

.. compounds-end

Reactions
.........

.. reactions-start


General format:
***************

.. code:: none 

   # Line that starts with "#" is considered as comment and ignored
   reactants    symbol    products    !energy    unit    scan
   reactants2   symbol    products2   !TSLabel
   reactants3   symbol    products3   !TSLabel
   TSLabel    energy    unit    scan

The first three examples after the comment correspond to the format of lines 
that indicate a reaction. The last one corresponds to the specification of a 
transition state, that can be either shared bewteen different reactions or 
just specified at the end for simplicity. 

Symbols:
********


+---------+--------------------------------------------------+
| Symbol  |   Description                                    | 
+=========+==================================================+
| **=>**  |   Forward reaction                               |
+---------+--------------------------------------------------+
| **<=**  |   Reverse reaction (When the --relative          |
|         |   flag is specified it uses as reference         |
|         |   the left side).                                |
+---------+--------------------------------------------------+
| **<=>** |   Forward and reverse reaction sharing           |
|         |   the same TS.                                   |
+---------+--------------------------------------------------+
| **<d>** |   Forward and reverse "diffusion" reaction.      |
|         |   Both share the TS. The specified energy is     |
|         |   used as: G(specified) = G(TS)                  |
|         |   - max(G(reactants),G(products))                |
+---------+--------------------------------------------------+


All symbols must be surrounded by empty spaces or tabs and should be used only 
once per reaction. 


Reversible reactions are automatically converted into 2 reactions,
the forward and the reverse reactions **sharing the same TS**. This is only 
important when setting up scans since: 

.. code:: none

   A  => B  !2.0 scan 
   B  => A  !2.0 

Does not behave as:

.. code:: none

   A <=> B !2.0 scan

However, the following reactions do have the same behaviour: 

.. code:: none 

   A  => B  !TS1
   B  => A  !TS1
   TS1  2.0  scan



Reactions file example: 
***********************

.. code:: none

   # Adduct Formation
    A         +    B     <d>    [A···B]                !2.0 kcal/mol
   [A···B]    +    Cat   <d>    [A···B···Cat]          !2.0 kcal/mol
   
   # Uncatalyzed reaction upon collision
    A         +    B     <=>     C                     !TS1
   
   # Uncatalyzed reaction with preformation of adduct
   [A···B]               <=>     C                     !TS1
   
   # Catalyzed reaction
   [A···B···Cat]         <=>    [C···Cat]              !10.0 kcal/mol scan
   
   # Product Diffusion
   [C···Cat]             <d>     C      +     Cat      !2.0 kcal/mol
   
   TS1     20.0 kcal/mol    scan

.. reactions-example-end

.. note::
   
   The --relative of pykinetic-model.py and pykinetic-scan.py will assume 
   that the energies specified are relative to the energy of the reactants. In 
   this example, enabling the --relative flag will result in the catalyzed 
   reaction having direct barrier of 10 kcal/mol (assuming a scan value of 0 
   kcal/mol). Not enabling the --relative flag will result in a barrier of 
   10.0 - E([A···B···Cat]) kcal/mol, with the Compounds file example it would 
   result in a barrier of 6.0 kcal/mol.

.. reactions-end

Simulation parameters
.....................

.. parameters-start

This file is optional and is generally useful when you already have a model
and are going to generate different versions of the model and you don not want 
to keep the tedious task of re-writing these parameters onto the generated file. 

Example file: 

.. code::  none

   dt          1E-12
   tfin        1E+02
   trep        0.1
   concentrations    0,1.0 ; 1,1.0 ; 2,0.1


*  The order may be changed.
*  Not all the parameters need to be specified in this file. 
*  dt is the solver timestep for calculating the solution. Some solvers require 
   it and others no.
*  trep  corresponds to the time step for reporting the concentrations of 
   all the species.  
*  tfin is the final time. Both number notations, 1E+02 and 100, are valid.
*  dt, trep and tfin are specified in seconds.
*  In the concentrations keyword only the non-zero initial concentrations in 
   M need to be included. In this example, the compounds 0 and 1 start with 1M
   and compound 2 with 0.1 M. (The number that corresponds to each compound 
   corresponds to the order in which they appear at the compounds file, 
   starting from 0)

.. parameters-end

Packaged Scripts Examples
-------------------------

pykinetic-model.py
..................

.. model-start

Creates a {python|c++} script that contains a system of mass balance equations 
(that represent a certain Chemical System) and the tools to solve it numerically.


.. code:: shell-session

   $ pykinetic-model.py compounds.txt reactions.txt model.py --writer python
   $ # After editing the file and including the concentrations and simulation times
   $ python model.py
   $ # A fast approach to visualize the data generated is with xmgrace
   $ # model.data is just a csv file where the first column is the time.
   $ xmgrace -nxy model.data

.. code:: shell-session

   $ pykinetic-model.py compounds.txt reactions.txt model.cpp --writer c++
   $ # After editing the file and including the concentrations and simulation times
   $ g++ model.cpp -o model.exe
   $ ./model.exe

.. model-end

pykinetic-scan.py
.................

.. scan-start

Creates a {python|c++} script that contains a system of mass balance equations 
(that represent a certain Chemical System) and the tools to solve it numerically
per each value of a scanned energy/energies of compounds/TSs. Currently it is 
recommended to do a --dryrun and run each generated model by itself.

.. code:: shell-session

   $ # it is recommended to use a simulation parameters file
   $ pykinetic-scan.py compounds.txt reactions.txt 0.0 0.5 5 --writer python --dryrun --I --simulation sim_params.txt


.. scan-end

Library usage examples
----------------------

.. library-usage-start

Creating a model from scratch 
.............................

Lets assume the following compounds and reactions: 

.. code:: none

   A    0.0  kcal/mol
   B    0.0  kcal/mol
   C    2.0  kcal/mol
   D   -2.0  kcal/mol

.. code:: none

   A + B <=> C   !10.0 kcal/mol
   C      => D   !18.0 kcal/mol

Now we create the model from scratch assuming a T of 25ºC. 

.. code:: python

   from pykinetic.classes import (ChemicalSystem, Energy, Compound, Reaction,
                                  TransitionState)
   from pykinetic.writers.python import Batch as Writer
   
   # We initialize the ChemicalSystem
   chemsys = ChemicalSystem(T=298.15)

   # Now we create the compounds
   A = Compound(A,Energy(0.0,'kcal/mol')) 
   B = Compound(B,Energy(1.0,'kcal/mol'))
   C = Compound(C,Energy(2.0,'kcal/mol'))
   D = Compound(D,Energy(-2.0,'kcal/mol'))
   
   # we can now add them one by one to the system:
   # chemsys.cadd(A)
   # chemsys.cadd(B) 
   # ...
   # or create a list and add them all
   compounds = [A,B,C,D] 
   chemsys.cextend(compounds,update=True)
   # if we dont update=True we will have to do it "manually"
   # chemsys.cupdate()

   # Now we create the reactions
   r1d = Reaction(reactants=(A,B),products=(C,)) 
   r1i = Reaction(reactants=(C,),products=(A,B))
   r2 = Reaction(reactants=(C,),products=(D,)) 
   # we can ignore the T parameter as it will be automatically set up by the 
   # ChemicalSystem object

   # We create the TSs and assign them to their reactions
   TS1 = TransitionState(Energy(10.0,'kcal/mol'),label='TS1',reactions=[r1d,r1i])
   TS2 = TransitionState(Energy(18.0,'kcal/mol'),label='TS2')
   TS2.reactions.append(r2)

   # Now we add the reactions to the ChemicalSystem. Now we use radd or rextend. 
   chemsys.radd(r1d)
   chemsys.radd(r1i)
   chemsys.radd(r2)

   # Now we have our chemical system already set up. Now we proceed to write it. 
   writer = Writer()
   writer.write(chemsys,'model.py')


Modifying and loading a model
.............................

.. code:: python

   from pykinetic.classes import (Energy, Compound, Reaction,
                                  TransitionState)
   from pykinetic.utils import BiasedChemicalSystem, calc_standard_state_correction
   from pykinetic.writers.python import Batch as Writer
   from pykinetic.userinput import populate_chemicalsystem_fromfiles
   
   # We initialize the ChemicalSystem and we want to apply a SS correction
   # from 1 atm -> 1 M
   std_correction = calc_standard_state_correction(T=298.15)
   chemsys = BiasedChemicalSystem(unit='kcal/mol',T=298.15)

   file_c = 'compounds.txt'
   file_r = 'reactions.txt'

   # Now we add from the files the compounds, reactions and TSs. 
   populate_chemicalsystem_fromfiles(chemsys,file_c,file_r,
                                      energy_unit='kcal/mol',
                                      relativeE=True)
   
   # Now we apply the bias. The biased system adds the bias to all the 
   # compounds and TSs. In this case it is applying the standard state correction
   # from 1 atm -> 1 M. 
   chemsys.apply_bias() 
   # The default bias is 0.0, but it is important to run this method after 
   # adding all the compounds and reactions when the bias is not 0.
   
   # we can write the model without the bias now 
   writer = Writer()
   writer.write(chemsys,'model.py')

   # Now we can change the bias if we decide so

   chemsys.bias = calc_standard_state_correction(T=298.15)

   # Now we proceed to write the model with std state correction. 
   writer.write(chemsys,'model_stdcorr.py')

.. library-usage-end