import argparse

# Initialize parser
# parser = argparse.ArgumentParser()
parser=argparse.ArgumentParser(
    description='''You did not specify any parameters. You must at least specify the number of chromosomes sampled and the sequence length. ''',
    epilog="""All's well that ends well.""")

# Adding optional argument

parser.add_argument('-n', "--TaxaNumber", type=int, default=10, help='Sets the number of isolates (default is 10)')
parser.add_argument('-g', "--Genome", type=int, default=5000, help='Sets the number and lengths of fragments of genetic material (default is 5000)')
parser.add_argument('-l', "--RecomLen", type=int, default=500, help='Sets the average length of an external recombinant interval, (default is 500)')
parser.add_argument('-r', "--RecomRate",type=float, default=0.05, help='Sets the site-specific rate of external (between species) recombination, (default is 0.05)')
parser.add_argument('-nu',"--nu" ,  type=float, default=0.2, help='nu')

# Read arguments from command line
args = parser.parse_args()

if args.TaxaNumber:
    print("Diplaying TaxaNumber as: % s" % args.TaxaNumber)

if args.Genome:
    print("Diplaying Alignment length as: % s" % args.Genome)

if args.RecomLen:
    print("Diplaying Recombination length as: % s" % args.RecomLen)

if args.RecomRate:
    print("Diplaying Recombination Rate as: % s" % args.RecomRate)

if args.nu:
    print("Diplaying nu as: % s" % args.nu)