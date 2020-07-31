#!/usr/bin/env python

import re
import sys
from argparse import ArgumentParser


def readParsnp11(input, cutoff, save_blocks):
    seqCountRE = re.compile(r'#SequenceCount\s+(\d+)')
    intervalCountRE = re.compile(r'#IntervalCount\s+(\d+)')
    seqIndexRE = re.compile(r'##SequenceIndex\s+(\d+)')
    seqFileRE = re.compile(r'##SequenceFile\s+(.+)')
    seqHeaderRE = re.compile(r'##SequenceHeader\s+(.+)')
    headerRE = re.compile(r'>(\d+):\d+-\d+\s+[\+-]\s+(\w+)')

    seqIndex = 0

    header = ''
    clusterID = ''

    headers = {}
    sequences = []
    concat = {}

    with open(input) as fp:
        for line in fp:
            line = line.rstrip('\n').rstrip('\r')
            if line.startswith('#'):
                seqCountMatch = seqCountRE.match(line)
                seqIntevalCountMatch = intervalCountRE.match(line)
                if seqCountMatch:
                    sys.stdout.write('Number of sequences {0}\n'.format(seqCountMatch.group(1)))
                if seqIntevalCountMatch:
                    sys.stdout.write('Number of blocks {0}\n'.format(seqIntevalCountMatch.group(1)))
                elif line.startswith('##'):
                    seqIndexMatch = seqIndexRE.match(line)
                    seqFileMatch = seqFileRE.match(line)
                    seqHeaderMatch = seqHeaderRE.match(line)
                    if seqIndexMatch:
                        seqIndex = seqIndexMatch.group(1)
                        headers[seqIndex] = {}
                    elif seqFileMatch:
                        headers[seqIndex]['FILE'] = seqFileMatch.group(1)
                    elif seqHeaderMatch:
                        headers[seqIndex]['HEADER'] = seqHeaderMatch.group(1).lstrip('>')
            else:
                if line.startswith('>'):
                    headerMatch = headerRE.match(line)
                    if headerMatch:
                        header = headerMatch.group(1)
                        clusterID = headerMatch.group(2)
                        sequences.append('')
                elif line.startswith('='):
                    if len(sequences[0]) >= cutoff:
                        print(clusterID)
                        if save_blocks:
                            with open(clusterID + '.fasta', 'w') as fasta_file:
                                for i in range(len(sequences)):
                                    fasta_file.write(
                                        '>' + headers[str(i + 1)]['FILE'] + '\n' + sequences[i].upper() + '\n')

                        if len(concat) == 0:
                            for i in range(len(sequences)):
                                concat[headers[str(i + 1)]['FILE']] = sequences[i].upper()
                        else:
                            for i in range(len(sequences)):
                                concat[headers[str(i + 1)]['FILE']] += sequences[i].upper()

                    sequences = []
                else:
                    sequences[int(header) - 1] += line
    return concat


def readMauve(input, cutoff, save_blocks):
    headerRE = re.compile(r'>\s*(\d+):\d+-\d+\s+[\+-]\s+(.+)')
    seqFileRE = re.compile(r'#Sequence(\d+)File\s+(.+)')
    headers = {}
    concat = {}
    sequences = {}
    counter = 0

    with open(input) as fp:
        for line in fp:
            line = line.rstrip('\n').rstrip('\r')
            if line.startswith('#'):
                seqFileMatch = seqFileRE.match(line)
                if seqFileMatch:
                    seqIndex = seqFileMatch.group(1)
                    headers[seqIndex] = seqFileMatch.group(2)
            else:
                if line.startswith('>'):
                    headerMatch = headerRE.match(line)
                    if headerMatch:
                        header = headerMatch.group(1)
                        sequences[header] = ''
                    else:
                        sys.exit('Could not read sequence #SequenceFile\n' + line)
                elif line.startswith('='):
                    aln_length = len(sequences.itervalues().next())
                    if aln_length >= cutoff:
                        if save_blocks:
                            with open('block{}.fna'.format(counter), 'w') as fasta_file:
                                for idx, seq in sequences.items():
                                    fasta_file.write('>' + headers[idx] + '\n' + seq.upper() + '\n')

                        for idx, seq in sequences.items():
                            concat[headers[idx]] = concat.get(headers[idx], '') + seq.upper()

                        for idx in set(headers.keys()) - set(sequences.keys()):
                            concat[headers[idx]] = concat.get(headers[idx], '') + ''.join(['-'] * aln_length);
                        counter += 1
                    sequences = {}
                else:
                    sequences[header] += line
    return concat


def main():
    parser = ArgumentParser(description='Convert an XMFA file to multiple FASTA files.')
    parser.add_argument('-i', '--input', metavar='FILE', required=True, help='an input sequence file in XMFA format.')
    parser.add_argument('-o', '--output', metavar='FILE', required=False,
                        help='an output sequence file in FASTA format containing concatenated sequences.')
    parser.add_argument('-c', '--cutoff', type=int, default=1, help='minimum alignment length.')
    parser.add_argument('-b', '--blocks', action="store_true", help='save each block to a separate fasta file')
    args = parser.parse_args()

    format = 'parsnp11'

    with open(args.input) as fp:
        for line in fp:
            if line.startswith('#FormatVersion Parsnp v1.1'):
                format = 'parsnp11'
                break
            elif line.startswith('#FormatVersion Mauve1'):
                format = 'mauve1'
                break

    if format == 'parsnp11':
        sequences = readParsnp11(args.input, args.cutoff, args.blocks)
    elif format == 'mauve1':
        sequences = readMauve(args.input, args.cutoff, args.blocks)
    else:
        sys.exit('File format not recognized')

    if args.output:
        with open(args.output, 'w') as fp:
            for h in sequences.keys():
                fp.write('>' + h + '\n' + sequences[h].upper() + '\n')


if __name__ == '__main__':
    main()