import argparse


parser=argparse.ArgumentParser(
    description='''You did not specify any parameters. You must at least specify the number of chromosomes sampled and the sequence length. ''',
    epilog="""All's well that ends well.""")
parser.add_argument('-n', type=int, default=10, help='Sets the number of isolates (default is 10)')
parser.add_argument('-g', type=int, default=5000, help='Sets the number and lengths of fragments of genetic material (default is 5000)')
parser.add_argument('-l', type=int, default=500, help='Sets the average length of an external recombinant interval, (default is 500)')
parser.add_argument('-r', type=float, default=0.05, help='Sets the site-specific rate of external (between species) recombination, (default is 0.05)')
parser.add_argument('-nu', type=float, default=0.2, help='nu')
args=parser.parse_args()