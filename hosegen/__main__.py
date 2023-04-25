from rdkit import Chem
import sys
from hosegen import HoseGenerator
from hosegen.geometry import *
import argparse

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

parser = argparse.ArgumentParser()
parser.add_argument("filename", help = "The molfile to use")
parser.add_argument("-a", "--AtomNumber", type=int, help = "Atom number to generate HOSE code for, starts at 0, if missing, all atoms will be done", required=False)
parser.add_argument("-s", "--Strict", type=str2bool, nargs='?', const=True, default=False, help = "strict mode (default: %(default)s)", required=False)
parser.add_argument("-u", "--UseStereo", type=str2bool, nargs='?', const=True, default=True, help = "generate stereo HOSE code (default: %(default)s)", required=False)
parser.add_argument("-m", "--MaxRadius", type=int, default=5, help = "number of spheres to use (default: %(default)s)", required=False)
parser.add_argument("-r", "--RingSize", type=str2bool, nargs='?', const=True, default=False, help = "include ring size in centre code (default: %(default)s)", required=False)
 
# Read arguments from command line
args = parser.parse_args()
 
path = args.filename
gen = HoseGenerator()
molStereo1 = Chem.MolFromMolFile(path, removeHs=False)
wedgemap1=create_wedgemap(path)
if args.AtomNumber is not None:
    value = gen.get_Hose_codes(molStereo1,args.AtomNumber,usestereo=args.UseStereo,strict=args.Strict,max_radius=args.MaxRadius,ringsize=args.RingSize,wedgebond=wedgemap1)
    print(value)
else:
    for idx,atom in enumerate(molStereo1.GetAtoms()):
        value = gen.get_Hose_codes(molStereo1,idx,usestereo=args.UseStereo,strict=args.Strict,max_radius=args.MaxRadius,ringsize=args.RingSize,wedgebond=wedgemap1)
        print(idx,":",value)

