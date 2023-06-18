=======
Writers
=======

Currently the main supported language is Python, but writers in C++ are also 
available. 

Python
------

Currently we have 4 possible reactor models for the Python Writer: 
Batch, SemiBatch, Extended SemiBatch, and Plug Flow Reactor (PFR). 

The simulation of the first three is with respect of time whereas the PFR's 
simulation is with respect to space. The SemiBatch reactor is mainly considered 
for experiments where a reactant is added dropwise, and its extended version 
accounts for a reactor that during some initial time has a reactant added 
dropwise and after the reactant has been exhausted it behaves as a Batch reactor.

C++
---

Currently only a single reactor model is available. The solution of kinetic 
models using this writer is discouraged, and it is aimed as an example of 
how to extend the code to other programing languages. 