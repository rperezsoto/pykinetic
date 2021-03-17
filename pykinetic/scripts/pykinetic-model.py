#!/usr/bin/python3
import os
import math
import argparse

from pathlib import Path

from pykinetic.classes import Energy,SimulationParameters,ConvergenceParameters
from pykinetic.writers import PythonWriter,CplusplusWriter
from pykinetic.utils import (BiasedChemicalSystem, write_indexfile, 
                              calc_standard_state_correction)
from pykinetic.userinput import populate_chemicalsystem_fromfiles

__version__ = "0.0.0"

WRITERS = {'python':PythonWriter,
           'c++':CplusplusWriter}

def create_parser():
    """ Create a Command line argument parser """
    description =""" Creates a python script that contains a
    system of differential equations that represent a Chemical System
    and the tools to solve it numerically"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("compounds",type=Path,
                        help="File with the compounds and Energies")
    parser.add_argument("reactions",type=Path, 
                        help="File with the reactions, energies and/or TSs")
    parser.add_argument("outfile", type=Path, help="Output script")
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
                        compounds and TSs""",
                        )
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

    return args


def main():
    parser = create_parser()
    args = parse_arguments(parser)
    unit = args.unit
    T = args.Temperature
    # Initialize the Writer
    writer_cls = WRITERS[args.writer]
    writer = writer_cls(conc_var='x',mb_var='dxdt',fun_var='model',
                        jac_var='Jac',jac_fun_var='jacobian',
                        header=args.header,tail=args.tail)
    writer.set_parameters(simulation=args.simulation,
                          convergence=args.convergence)
    # Initialize the ChemicalSystem 
    chemsys = BiasedChemicalSystem(bias=args.bias,T=args.Temperature,unit=unit)
    populate_chemicalsystem_fromfiles(chemsys,
                                      file_c=args.compounds,
                                      file_r=args.reactions,
                                      energy_unit=unit,
                                      relativeE=args.relative)
    chemsys.apply_bias()

    data_file_name = f'{args.outfile.stem}.data' 
    writer.parameters['out_filename'] = data_file_name
    writer.write(chemsys,args.outfile)

    if args.IndexFile:
        IdxFile = args.outfile.parent.joinpath(args.outfile.stem + '.index')
        write_indexfile(chemsys,IdxFile,isrelative=args.relative)

if __name__ == '__main__':
    main()
