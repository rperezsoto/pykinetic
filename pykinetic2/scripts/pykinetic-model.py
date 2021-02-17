import sys
import os
import math
import argparse

from pykinetic2.Classes import ChemicalSystem, Energy
from pykinetic2.Writers import PythonWriter,CplusplusWriter
from pykinetic2.Utils import BiasedChemicalSystem, SimulationParameters, ConvergenceParameters
from pykinetic2.InputParse import populate_chemicalsystem_fromfiles

__version__ = "0.0.0"

def CreateParser():
    """ Create a Command line argument parser """
    Description =""" Creates a python script that contains a
    system of differential equations that represent a Chemical System
    and the tools to solve it numerically"""
    parser = argparse.ArgumentParser(description=Description)
    parser.add_argument("compounds",help="File with the compounds and Energies")
    parser.add_argument("reactions",help="File with the reactions and Energies")
    parser.add_argument("OFile",help="File with the reactions and Energies")
    parser.add_argument("-E","--EParse",choices=["absolute","relative"],
                type=lambda s:s.lower(),default='absolute',
                help="""Whether the energy in the reactions file is relative to
                reactants or the absolute value of the TS structure""")
    parser.add_argument("-I","--IndexFile",help="""If Enabled creates a txt file
                        named 'OutFile'.index with the indices used for
                        reactions and compounds """,
                        default=False,action="store_true")
    parser.add_argument("--dot",help="""If Enabled creates a .dot file named
                        ready for its coversion to an image with the 'dot' tool""",
                        default=False,action="store_true")
    parser.add_argument("-T","--Temperature", nargs=2,metavar=("Temp","Unit"),
                        default=(25.0,'C'),help="""Temperature and unit
                        ('C'=Celsius and 'K'=Kelvin) """)
    parser.add_argument("--Unit",help="""Energy unit used in the Input Files
                        Diffusion reactions '<d>' are treated differently and always
                        assumed to be given in kcal/mol""", default='kcal/mol',
                        choices='hartree eV cm-1 kcal/mol kJ/mol J/mol K J Hz'.split())
    parser.add_argument("--Correction",help="""Applies a correction to all energies.
                        if no argument is provided: applies Standard State correction
                        otherwise uses the specified value (kcal/mol) """,
                        nargs='?',default=0.0, const='Standard State')
    parser.add_argument("--Conv_Params", help="""Path to the file with the
                        Convergence parameters""", default=None)
    parser.add_argument("--Sim_Params", help="""Path to the file with the
                        Simulation parameters""", default=None)
    parser.add_argument("--Tail",help="""Path to the non-standard Template for
                        the ending part of the generated script""",default=None)
    parser.add_argument("--Header",help="""Path to the non-standard Template for
                        the head of the generated script""",default=None)
    parser.add_argument("--Function",help="""Path to the non-standard function
                        structure for the generated script""",default=None)
    parser.add_argument("--Jacobian",help="""Path to the non-standard Jacobian
                        structure for the generated script""",default=None)
    return parser
def ParseArguments(parser):
    """ Check the arguments of the parser and prepare them for their use """
    args = parser.parse_args()
    args.compounds = os.path.abspath(args.compounds)
    args.reactions = os.path.abspath(args.reactions)
    args.OFile = os.path.abspath(args.OFile)
    if args.Conv_Params is not None:
        args.Conv_Params = os.path.abspath(args.Conv_Params)
    if args.Sim_Params is not None:
        args.Sim_Params = os.path.abspath(args.Sim_Params)
    if args.Header is not None:
        args.Header = os.path.abspath(args.Header)
    if args.Tail is not None:
        args.Tail = os.path.abspath(args.Tail)
    if args.Function is not None:
        args.Function = os.path.abspath(args.Function)
    if args.Jacobian is not None:
        args.Jacobian = os.path.abspath(args.Jacobian)
    if args.Temperature[1] == 'C':
        T = float(args.Temperature[0]) + 273.15
    elif args.Temperature[1] == 'K':
        T = float(args.Temperature[0])
    else:
        msg = "Temperature in {} not implemented".format(args.Temperature[1])
        raise NotImplementedError(msg)
    args.Temperature = T
    if args.Correction == 'Standard State':
        args.Correction = Energy(8.314*T*math.log(0.082*T),'J/mol')
    else:
        args.Correction = Energy(args.Correction,'kcal/mol')
    return args

def main():
    parser = CreateParser()
    args = ParseArguments(parser)
    # Initialize the ChemicalSystem
    T, bias, unit = args.Temperature, args.Correction, args.Unit
    writer = CplusplusWriter(conc_var='x',mb_var='dxdt',fun_var='model',
                          jac_var='Jac',jac_fun_var='Jacobian',
                          header=args.Header,tail=args.Tail)
    sim_params = SimulationParameters.read_from(args.Sim_Params)
    if args.Conv_Params is None:
        conv_params = ConvergenceParameters()
    else:
        conv_params = ConvergenceParameters.read_from(args.Conv_Params)
    writer.set_parameters(simulation=sim_params,convergence=conv_params)
    chemsys = BiasedChemicalSystem(bias=bias,T=T,unit=unit)
    is_relativeE = args.EParse == 'relative'
    populate_chemicalsystem_fromfiles(chemsys,
                                      file_c=args.compounds,
                                      file_r=args.reactions,
                                      energy_unit=unit,
                                      relativeE=is_relativeE)
    writer.write(chemsys,args.OFile)

    if args.dot:
        Name = args.OFile[:-3]+'.dot'
        with open(Name,'w') as F :
            F.write('\n'.join(DotFile(ChemSys)))
    if args.IndexFile:
        IdxFile = args.OFile[:-3]+'.index'
        Write_IndexFile(ChemSys,IdxFile)

if __name__ == '__main__':
    main()
