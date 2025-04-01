import argparse
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--FilePath', type=str, default='1',
                        help='File path for MD trajectory')
    parser.add_argument('--Format', type=str, default='vasp-xdatcar',
                        help='File format for MD trajectory. Available formats: https://wiki.fysik.dtu.dk/ase/ase/io/io.html#ase.io.read')
    parser.add_argument('--StartFrame', type=int, default=0, 
                        help='Frames that should be deleted, index of start frame')
    parser.add_argument('--EndFrame', type=int,  default=-1,
                        help='Frames that should be deleted, index of end frame')
    parser.add_argument('--ExceptAtom', nargs='+', default='C',
                        help='Atoms that should be deleted')
    parser.add_argument('--PBC', nargs=3, default=[True,True,False], help='PBC for system')
    parser.add_argument('--InitialRadius', nargs='+', type=str, default='H=0.775 O=0.175',
                        help='The initial radius of a persistance homology, format of key-value pairs: key1=value1 key2=value2')
    return parser.parse_args()
def parse_key_value_pairs(items):
    result = {}
    items = items.split()
    for item in items:
        key, value = item.split("=")
        result[key] = float(value)
    return result

if __name__ == "__main__":
    args = parse_args()
    config_dict = parse_key_value_pairs(args.InitialRadius)
    print(config_dict)