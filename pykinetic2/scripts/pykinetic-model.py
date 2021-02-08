import sys
import os
import math
import argparse

from pykinetic.Classes import Compound, ElementalStep, Energy
from pykinetic.Utils import *

__version__ = "1.2.0"

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
    T, Corr, Unit = args.Temperature, args.Correction, args.Unit
    ChemSys = CorrChemSys(Corr,T,Unit)
    ChemSys.compounds_fromFile(args.compounds)
    ChemSys.reactions_fromFile(args.reactions,args.EParse)
    # Check and print for the user to correct the compounds.txt
    missing_compounds = ChemSys.check_compounds()
    if missing_compounds:
        print("The Following compounds are not defined in the compounds.txt")
        for item in sorted(missing_compounds):
            print('{}\tNUMBER'.format(item))
        raise SystemExit("Aborting program, Not all species properly specified")

    Template = FileTemplate(Sim_Params=args.Sim_Params,
                            Conv_Params=args.Conv_Params,
                            Header=args.Header,
                            Function=args.Function,
                            Jacobian=args.Jacobian,
                            Tail=args.Tail)
    Template.update_with(ChemSys)
    Template.Fill()
    Template.write_to(args.OFile)

    if args.dot:
        Name = args.OFile[:-3]+'.dot'
        with open(Name,'w') as F :
            F.write('\n'.join(DotFile(ChemSys)))
    if args.IndexFile:
        IdxFile = args.OFile[:-3]+'.index'
        Write_IndexFile(ChemSys,IdxFile)

if __name__ == '__main__':
    main()
