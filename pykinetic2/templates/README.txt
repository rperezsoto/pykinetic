These are the templates for the writing of the Script that can perform the
simulation.

Some rules for writing templates have to be taken into account:

-A workaround due to the .format python string setting, to write a "{}" in the
output file write instead "{0.fill}"

-A workaround to the reading of template files to include a final \n is to
write a "{0.end}"

- Be careful and ensure that the template files have \t instead of leading \b
although it shouldn't give problems, it's better to be careful.  (due to the
line.strip() routine used when reading)

- initial concentrations of each compound are weird to write. As a general
rule one won't have a very large number of these so writing
>> Conc_ini	1,1.0 ; 2,0.5 ; 34,2.0 ; [...]
although not the most beautiful thing, is the simplest to implement without
having to split the file Sym_Parameters into two.
** And please, spaces around ; are OK but ensure that the i,j are not
separated with spaces. Will not change the result nor break the code.
** Also please, for i,j i=integer and j=float. be explicit with this, not
respecting this can actually change the result of the simulation or raise an
error when trying to execute it.

Bests,

The developer
