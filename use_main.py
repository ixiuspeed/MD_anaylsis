from use_argparse import parse_args
from ase.io import read,write
from MD_anaylsis import Anaylsis



if __name__ == '__main__':
    args = parse_args()
    frame = read(filename=args.FilePath, index=':',format=args.Format)
    frame = frame[args.StartFrame:args.EndFrame]
    for atoms in frame:
        del atoms[[atom.index for atom in atoms if (atom.symbol in args.ExceptAtom)]]
        atoms.set_pbc(args.PBC)
        anayls = Anaylsis(atoms, args)
