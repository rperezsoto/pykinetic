#!/usr/bin/python
import os
import math
import argparse
import subprocess

from pathlib import Path

from pykinetic2.Classes import Energy,SimulationParameters,ConvergenceParameters
from pykinetic2.Writers import PythonWriter,CplusplusWriter
from pykinetic2.Utils import ScannableChemicalSystem,write_indexfile
from pykinetic2.InputParse import populate_chemicalsystem_fromfiles

def CreateParser():
    """ Create a Command line argument parser """
    description =""" Creates a python script that contains a
    system of differential equations that represent a Chemical System
    and the tools to solve it numerically"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("compounds",type=Path,
                        help="File with the compounds and Energies")
    parser.add_argument("reactions",type=Path, 
                        help="File with the reactions, energies and/or TSs")
    parser.add_argument("start",help="Starting value for the correction scan")
    parser.add_argument("stop",help="Final value for the correction scan")
    parser.add_argument("steps",help="Number of steps in between")
    parser.add_argument("--scanunit",help="Unit of Start and Stop parameters",
                        default='kcal/mol',
                        choices=Energy._units)
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
    parser.add_argument("--ScriptFiles",help="""If Enabled retains the generated
                        .py files""",default=False,action="store_true")
    parser.add_argument("--dryrun",default=False,action="store_true",
                        help="""If enabled it will generate the scripts but wont
                        run them""")
                        
    return parser
def ParseArguments(parser):
    """ Check the arguments of the parser and prepare them for their use """
    args = parser.parse_args()
    args.compounds = os.path.abspath(args.compounds)
    args.reactions = os.path.abspath(args.reactions)
    args.Start = float(args.Start)
    args.Stop = float(args.Stop)
    args.Steps = int(args.Steps)
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
    if args.Std_State:
        args.Std_State = Energy(8.314*T*math.log(0.082*T),'J/mol')
    if args.dryrun:
        args.ScriptFiles = True
    return args


def main():
    parser = CreateParser()
    args = ParseArguments(parser)
    # Initialize the ChemicalSystem
    T, Std_State, Unit = args.Temperature, args.Std_State, args.Unit
    SUnit, Start, Stop, Steps = args.SUnit, args.Start, args.Stop, args.Steps
    step = (Stop - Start)/Steps
    Steps += 1 # in order to include the final point.
    E_step = Energy(step,SUnit)
    E_ini = Energy(Start,SUnit)
    tmp_file = 'Temp.py'
    if not args.dryrun:
        mkdir_p('Outputs')
    if args.IndexFile:
        mkdir_p('IndexFiles')
    if args.ScriptFiles:
        mkdir_p('ScriptFiles')
    if Std_State:
        E_ini = E_ini + Std_State
    # Record in a file the value of the correction of each step
    int_f = '{:03d}'.format
    energy_f = lambda x: str(x).split()[0]
    with open('Scan_Params.out.txt','w') as F:
        F.write('Scan Step\tEnergy({})\n'.format(args.SUnit))
    # Create ChemSys
    ChemSys = MutableSystem(E_ini,T,Unit)
    ChemSys.compounds_fromFile(args.compounds)
    ChemSys.reactions_fromFile(args.reactions,args.EParse)
    Template = FileTemplate(Sim_Params=args.Sim_Params,
                            Conv_Params=args.Conv_Params,
                            Header=args.Header,
                            Function=args.Function,
                            Jacobian=args.Jacobian,
                            Tail=args.Tail)
    for i in range(Steps):
        Template.update_with(ChemSys)
        # Create Name of Script and OutFile
        Name = 'Outputs/Scan_{:03d}.out'.format(i)
        print(Name)
        Template.OFileName = Name
        Template.Fill()
        Template.write_to(tmp_file)
        if not args.dryrun:
            # Run Script as subprocess
            cmd = ['python', tmp_file]
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            #for line in p.stdout:
            #    print(line)
            # Wait for it to finish
            p.wait()
        if args.IndexFile:
            IdxFile = 'IndexFiles/Scan_{:03d}.index'.format(i)
            Write_IndexFile(ChemSys,IdxFile)
        if args.ScriptFiles:
            os.rename(tmp_file,'ScriptFiles/Scan_{:03d}.py'.format(i))
        with open('Scan_Params.out.txt','a') as F:
            F.write('{: ^9}\t{: ^16}\n'.format(int_f(i),energy_f(E_step*i + Start)))
        # Update Energy
        ChemSys.MutateEnergy(E_step)

if __name__ == "__main__":
    main()
