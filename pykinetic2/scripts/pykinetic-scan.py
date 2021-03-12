#!/usr/bin/python3
import os
import math
import argparse
import subprocess

from pathlib import Path

from pykinetic.classes import Energy,SimulationParameters,ConvergenceParameters
from pykinetic.writers import PythonWriter,CplusplusWriter
from pykinetic.utils import (ScannableChemicalSystem, write_indexfile,
                              calc_standard_state_correction)
from pykinetic.userinput import populate_chemicalsystem_fromfiles

__version__ = "0.0.0"

WRITERS = {'python':PythonWriter,
           'c++':CplusplusWriter}
COMMANDS = {}
COMMANDS['python'] = lambda x: ['python3',x ]
COMMANDS['c++'] = lambda x: ['g++',str(x),';',f'./{x}']

def create_parser():
    """ Create a Command line argument parser """
    description =""" Creates a python script that contains a
    system of differential equations that represent a Chemical System
    and the tools to solve it numerically"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("compounds",
                        type=Path,
                        help="File with the compounds and Energies")
    parser.add_argument("reactions",
                        type=Path, 
                        help="File with the reactions, energies and/or TSs")
    parser.add_argument("start",
                        type=float,
                        help="Starting value for the correction scan")
    parser.add_argument("stop",
                        type=float,
                        help="Final value for the correction scan")
    parser.add_argument("steps",type=int,
                        help="Number of steps in between")
    parser.add_argument("--scanunit",
                        default='kcal/mol',
                        choices=Energy._units,
                        help="Unit of Start and Stop parameters")
    parser.add_argument("outfile", # No deberia existir a priori
                        type=Path, 
                        help="Output script")
    parser.add_argument("--relative",
                        action='store_true',
                        default=False,
                        help="""Whether the energy in the reactions file is 
                        relative to reactants. Otherwise the absolute value of 
                        the TS structure is assumed to correspond to the 
                        provided energy""")
    parser.add_argument("-I","--IndexFile",
                        action="store_true",
                        default=False,
                        help="""If Enabled creates a txt file named 
                        'OutFile'.index with the indices used for reactions,
                        compounds and TSs""")
    parser.add_argument("-T","--Temperature", 
                        nargs=2,
                        metavar=("Temp","Unit"),
                        default=(25.0,'C'),
                        help="Temperature and unit ('C'=Celsius and 'K'=Kelvin)")
    parser.add_argument("--unit",
                        default='kcal/mol',
                        help="""Energy unit used in all the energy values whose 
                        unit is not specified""", 
                        choices=Energy._units)
    parser.add_argument("--bias",
                        nargs='?',
                        default=0.0,
                        const='Standard State',
                        help="""Applies a bias to all energies. If no argument 
                        is provided it applies Standard State correction. 
                        Otherwise, uses the specified value """)
    parser.add_argument("--convergence", 
                        type=Path,
                        default=None,
                        help="Path to the file with the Convergence parameters")
    parser.add_argument("--simulation",
                        type=Path,
                        default=None,
                        help="Path to the file with the Simulation parameters")
    parser.add_argument("--tail",
                        type=Path,
                        default=None,
                        help="""Path to the non-standard Template for the ending
                         part of the generated script""")
    parser.add_argument("--header",
                        type=Path,
                        default=None,
                        help=""" Path to the non-standard Template for the head 
                        of the generated script""")
    parser.add_argument("--writer",
                        choices=['python','c++'],
                        default='python',
                        help="format for the output model")
    parser.add_argument("--scripts",
                        default=False,
                        action="store_true",
                        help="""If Enabled retains the generated scripts""")
    parser.add_argument("--dryrun",
                        default=False,
                        action="store_true",
                        help="""If enabled it will generate the scripts but wont
                        run them""")
                        
    return parser
def parse_arguments(parser):
    """ Check the arguments of the parser and prepare them for their use """
    args = parser.parse_args()

    # Read the header if it has been specified
    if args.header is not None: 
        with open(args.header) as F:
            txt = F.read()
        args.header = txt
    # Read the tail if it has been specified
    if args.tail is not None: 
        with open(args.tail) as F:
            txt = F.read()
        args.tail = txt

    # Handle Temperature
    if args.Temperature[1] == 'C':
        T = float(args.Temperature[0]) + 273.15
    elif args.Temperature[1] == 'K':
        T = float(args.Temperature[0])
    else:
        msg = "Temperature in {} not implemented".format(args.Temperature[1])
        raise NotImplementedError(msg)
    args.Temperature = T

    # Handle bias
    if args.bias == 'Standard State':
        bias = calc_standard_state_correction(T)
        bias.to_unit(args.unit)
    else:
        bias = Energy(args.bias,args.unit)
    args.bias = bias

    # Handle simulation parameters
    if args.simulation is None: 
        sim_params = SimulationParameters()
    else:
        sim_params = SimulationParameters.read_from(args.simulation)
    args.simulation = sim_params

    # Handle convergence parameters
    if args.convergence is None:
        conv_params = ConvergenceParameters()
    else:
        conv_params = ConvergenceParameters.read_from(args.convergence)
    args.convergence = conv_params

    if args.dryrun:
        args.scripts = True

    return args

def main():
    parser = create_parser()
    args = parse_arguments(parser)
    unit = args.unit
    scan_unit = args.scanunit
    # Prepare Filesystem
    outputs_folder = Path('Outputs')
    indices_folder = Path('IndexFiles')
    scripts_folder = Path('ScriptFiles')
    scan_parameters_file = Path('Scan_Params.out.txt')
    tmp_name = 'Temp'
    if args.writer == 'python':
        suffix = 'py'
    elif args.writer == 'c++':
        suffix = 'cpp'
    tmp_file = Path(f'{tmp_name}.{suffix}')
    if not args.dryrun:
        outputs_folder.mkdir(exist_ok=True)
    if args.IndexFile:
        indices_folder.mkdir(exist_ok=True)
    if args.scripts:
        scripts_folder.mkdir(exist_ok=True)

    # Initialize the Writer
    writer_cls = WRITERS[args.writer]
    writer = writer_cls(conc_var='x',mb_var='dxdt',fun_var='model',
                        jac_var='Jac',jac_fun_var='Jacobian',
                        header=args.header,tail=args.tail)
    writer.set_parameters(simulation=args.simulation,
                          convergence=args.convergence)

    # Initialize the ChemicalSystem
    step = (args.stop - args.start)/args.steps
    steps = args.steps + 1 # in order to include the final point.
    energies = [Energy(args.start + step*i,args.scanunit) for i in range(steps)]

    # Record in a file the value of the correction of each step
    int_f = '{:03d}'.format
    energy_f = lambda x: str(x).split()[0]
    with open(scan_parameters_file,'w') as F:
        F.write(f'Scan Step    Energy({args.scanunit})\n')
    
    # Initialize the ChemicalSystem 
    chemsys = ScannableChemicalSystem(scan_unit=scan_unit,unit=unit,
                                      bias=args.bias,T=args.Temperature)
    populate_chemicalsystem_fromfiles(chemsys,
                                      file_c=args.compounds,
                                      file_r=args.reactions,
                                      energy_unit=unit,
                                      relativeE=args.relative)
    chemsys.apply_bias()
    chemsys.apply_scan()

    for i,energy in enumerate(energies): # Keep going
        stem = f'scan_{i:03d}' # stem of the produced files
        chemsys.scan = energy

        # Create Name of Script and OutFile
        out_data_file = f'{stem}.data'
        writer.parameters['out_filename'] = out_data_file
        print(stem)

        # Write the script
        writer.write(chemsys,tmp_file)

        if not args.dryrun:
            # Run Script as subprocess
            cmd = COMMANDS[args.writer](tmp_file)
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            #for line in p.stdout:
            #    print(line)
            # Wait for it to finish
            p.wait()

        if args.IndexFile:
            index_file = indices_folder.joinpath(f'{stem}.index')
            write_indexfile(chemsys,index_file, isrelative=args.relative)
        
        if args.scripts:
            script = scripts_folder.joinpath(f'{stem}{tmp_file.suffix}')
            tmp_file.rename(script)

        with open(scan_parameters_file,'a') as F:
            F.write('{: ^9}\t{: ^16}\n'.format(f'{i:03d}',str(energy)))

if __name__ == "__main__":
    main()
