import argparse as ap
parser = ap.ArgumentParser(fromfile_prefix_chars="@")
parser.add_argument("--foo", default="AAA")
args=parser.parse_args()
print(args.foo)
