from compsci260lib import *


def run_assembly_tester():
    """Read in the read and contig files, verify the reads in each contig, and estimate
    the gap length between the contigs.
    """
    reads_file = "paired.reads.fasta"
    contigs_file = "contigs.fasta"

    # Load the two fasta files
    reads = get_fasta_dict(reads_file)
    contigs = get_fasta_dict(contigs_file)

    # Determine how the reads map to the contigs
    contig_reads = find_contig_reads(reads, contigs)

    # Report:
    #   - Any reads that are not consistent with the overall assembly based
    #     on the checks done by find_contig_reads()
    #   - The read pairs where both reads match to (distinct) contigs in supercontig1,
    #   and their estimated gap lengths
    #
    # report the non-consistent values by finding the ones that have None as the values of the dictionairy
    print("These are the reads not consistent with overall assembly:")
    valWNone = []
    contPairs = {}
    for cr in contig_reads:
        valList = list(contig_reads[cr].values())
        if None in valList:
            valWNone.append(cr)
        else:
            contPairs[cr] = ((valList[0], valList[3]))
    print(valWNone)
    print("There are a total of " + str(len(valWNone)) + " reads like this.")
    print("")

    # find the pairs of contigs with the found sequences that have contigs accross multiple values
    contPairsDiff = {}
    for x in contPairs:
        val = contPairs[x]
        if val[0] != val[1]:
            contPairsDiff[x] = (val)

    # find and report all the supercontigs
    print("There are " + str(len(contPairsDiff)) + " supercontigs")
    a = dict(sorted(contPairsDiff.items(), key=lambda item: item[1]))

    print("The supercontigs are (sorted by start contig):")
    for x in a:
        print(x, a[x])

    superContig = {}
    for x in a:
        if a[x][0] in superContig.keys():
            superContig[a[x][0]].append(x)
        else:
            superContig[a[x][0]] = []
            superContig[a[x][0]].append(x)

    print("")
    print("There are " + str(len(superContig)) + " supercontigs")
    print("The sequences in supercontig1 are: ")
    print(superContig['contig1'])

    # finding the gap sizes and create a table
    gaps = []
    for x in superContig['contig1']:
        length = 2000 - ((len(contigs['contig1'])-contig_reads[x]['start_a']) + contig_reads[x]['end_b'] +1)
        pm = str(length) + " Â± 10"
        gaps.append([x, pm, length-10, length+10])

    print("")
    print ("{:<10} {:<20} {:<15} {:<15}".format('Read Name', 'Gap Size Range', 'Lower Bound', 'Upper Bound'))
    for v in gaps:
        name, size, bottom, up = v
        print ("{:<10} {:<20} {:<15} {:<15}".format(name, size, bottom, up))

    #


def find_contig_reads(reads, contigs):
    """
    Determine the contig in which each sequencing read appears (if any), along with
    where in the contig it matches.  The `contigs` dict will have a collection of
    contig names mapped to contig sequences.  The `reads` dict contains separate keys
    (labeled 'a' and 'b') for the two reads in each read pair.  However, this function
    returns a dictionary with a single key for each read pair (i.e., without 'a' or 'b').

    The value for each read-pair key will itself be a dictionary with the following keys:

        - 'contig_a' (str): the contig in which read 'a' was found, as
          <contig name> or None (<contig name> will be a key in `contigs` dict)
        - 'start_a' (int): the start position (1-indexed) read 'a' mapped to
          within its respective contig (None if not found in any contig)
        - 'end_a' (int): the end position (1-indexed) read 'a' mapped to
          within its respective contig (None if not found in any contig)
        - 'contig_b' (str): the contig in which read 'b' was found, as
          <contig name> or None. (see 'contig_a' for example).
        - 'start_b' (int): the start position (1-indexed) read 'b' mapped to
          within its respective contig (None if not found in any contig)
        - 'end_b' (int): the end position (1-indexed) read 'b' mapped to
          within its respective contig (None if not found in any contig)

    The returned dictionary should look something like:
    {
        'seq1': {
            'contig_a': 'contig1',
            'start_a': 301,
            'end_a': 800,
            'contig_b': None
            'start_b': None,
            'end_b': None,
        },
        'seq2': {
            'contig_a': 'contig2',
            'start_a': 1101,
            'end_a': 1600,
            'contig_b': 'contig1'
            'start_b': 201,
            'end_b': 700,
        },
        'seq3' : {
            'contig_a': None,
            'start_a': None,
            'end_a': None,
            'contig_b': None
            'start_b': None,
            'end_b': None,
        },
        ...
    }

    Arguments:
        reads (dict str to str): dictionary of reads, mapping read names to sequences
        contigs (dict str to str): dictionary of contigs, mapping contig names to sequences

    Returns:
        Dictionary mapping read-pairs to information about their reads' locations in contigs.
    """
    retDict = {}


    for read in reads:
        # add the read to the dictionairy
        if read[:-1] not in retDict.keys():
            retDict[read[:-1]] = {}
        letter = read[len(read) - 1:] # get the letter for the specific sequence read
        # check if the read is in any of the contigs
        for contig in contigs:
            largeContig = contigs[contig]
            fRead = reads[read]
            f = largeContig.find(fRead) # find the full read sequence in the full contig sequence
            # if there is a match found, add the appropriate values to the dicitonairy
            if f != -1:
                retDict[read[:-1]]['contig_' + letter] = contig
                retDict[read[:-1]]['start_' + letter] = f + 1
                retDict[read[:-1]]['end_' + letter] = f + len(fRead)
        # if there is no values or matches found with the read for the contigs
        if ('contig_'+letter) not in retDict[read[:-1]]:
            retDict[read[:-1]]['contig_' + letter] = None
            retDict[read[:-1]]['start_' + letter] = None
            retDict[read[:-1]]['end_' + letter] = None

    return retDict


if __name__ == '__main__':
    """Call run_assembly_tester(). Do not modify this block."""
    run_assembly_tester()