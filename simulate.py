import sys
import random
from compsci260lib import *


def run_simulate():
    """
    Simulates the sequencing process then empirically compute and report
    the quantities from parts a-c.
    """

    iterations = 20
    G = 3000000
    R = 40000
    L = 450

    # Call simulate(G, R, L) `iteration` times. Report the empirical values
    # from parts a-c for each iteration
    #
    values = {'Empirical Coverage': [], 'Number of Nucleotides': [], 'Number of Contigs': [], 'Average Length': []}
    for x in range(iterations):
        print("Iteration " + str(x + 1))
        (EmpiricalCoverage, NumNucleotides, NumContigs, AvgLeng) = simulate(G, R, L)
        values['Empirical Coverage'].append(EmpiricalCoverage)
        values['Number of Nucleotides'].append(NumNucleotides)
        values['Number of Contigs'].append(NumContigs)
        values['Average Length'].append(AvgLeng)
        print(EmpiricalCoverage, NumNucleotides, NumContigs, AvgLeng)
        print("")

    # report the average values for all the iterations
    print("Average Values for iterations: ")
    print("Average Empirical Coverage: " + str(sum(values['Empirical Coverage']) / len(values['Empirical Coverage'])))
    print("Average Number of Nucleotides: " + str(
        sum(values['Number of Nucleotides']) / len(values['Number of Nucleotides'])))
    print("Average Number of Contigs: " + str(sum(values['Number of Contigs']) / len(values['Number of Contigs'])))
    print("Mean Average Length: " + str(sum(values['Average Length']) / len(values['Average Length'])))

    values.clear()
    #


def simulate(G, R, L):
    """
    Simulates one iteration of the sequencing process and empirically compute
    the empirical coverage (average number of times a nucleotide in the genome
    was sequenced), the number of nucleotides not covered by any read, the
    number of contigs assuming you can use the oracular assembly algorithm to
    assemble all the reads, and the average length of these contigs.

    Args:
        G (int) - the length of the genome
        R (int) - the number of reads
        L (int) - the length of each read

    Returns
        a tuple of floats:

            (Empirical coverage,
             Number of nucleotides not covered by any read,
             Number of contigs,
             Average length of these contigs)
    """

    #
    # initialize the list to 0 based on the size of the sequence
    listNucleotide = [[0] for i in range(G)]

    # get the start and fill in the list based on the times the nucleotide is sequenced
    for x in range(R):
        randomStart = random.randint(0, len(listNucleotide) - 1)
        for y in range(L):
            if (randomStart + y) < len(listNucleotide):
                listNucleotide[randomStart + y].append(listNucleotide[randomStart + y][0] + 1)
                listNucleotide[randomStart + y].pop(0)

    # calculate all the statistics by going through the list of sequenced nucleotides
    dictNumbers = {}
    dictNumbers['numberOfContigs'] = 0
    dictNumbers['length'] = []
    count = 0
    totalCount = 0
    zeroCount = 0

    # go through the sequence
    for x in range(len(listNucleotide)):
        if listNucleotide[x][0] != 0:
            if count == 0 and x != 0:
                dictNumbers['numberOfContigs'] += 1
            count += 1
            totalCount += listNucleotide[x][0]
        else:
            dictNumbers['length'].append(count)
            count = 0
            zeroCount += 1

    avg = sum(dictNumbers['length']) / dictNumbers['numberOfContigs']
    #
    return ((totalCount / len(listNucleotide)), zeroCount, dictNumbers['numberOfContigs'], avg)


if __name__ == '__main__':
    """Call run_simulate(), do not modify"""
    run_simulate()