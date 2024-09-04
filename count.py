import pandas as pd
import collections


def find_first_deletion(cigar):
    """Find the position of the the first deletion in a
    read from the cigar string, will return 0 if no deletion
    found"""

    position = 0
    for operation, length in cigar:
        if operation == 2:
            return position
        else:
            position += length

    return position


##################################################
def getCrosslink(read):
    """Finds the crosslinked base from a pysam read.

    Cross linked bases are definated as in Sugimoto et al, Genome Biology 2012

        The nucleotide preceding the iCLIP cDNAs mapped by Bowtie was used to
        define the cross link sites identified by truncated cDNAs.

        [For reads with deletions] The deleted nucleotide in CLIP and iCLIP
        cDNAs mapped by Novoalign was used to define the cross-link sites
        identified by read-through cDNAs. If a cDNA had more than one deletion,
        we selected the one closest to the beginning of the read.

    returns a tuple with the position of the read and one of the following
    categories:

        * truncated_neg

        * truncated_pos

        * deletion_neg

        * deletion_pos


    to record whether the position came from a truncation or a deletion"""

    if "D" not in read.cigarstring:
        if read.is_reverse:
            pos = read.aend

        else:
            pos = read.pos - 1

    else:
        if read.is_reverse:
            cigar = reversed(read.cigar)
            position = find_first_deletion(cigar)
            pos = read.aend - position - 1

        else:
            position = find_first_deletion(read.cigar)
            pos = read.pos + position

    return pos


##################################################
def countChr(reads, chr_len, dtype="uint16"):
    """Counts the crosslinked bases for each read in the pysam rowiterator
    reads and saves them in pandas Series: those on the positive strand
    and those on the negative strand. The Series are indexed on genome position,
    and are sparse.

    Cross linked bases are definated as in Sugimoto et al, Genome Biology 2012

        The nucleotide preceding the iCLIP cDNAs mapped by Bowtie was used to
        define the cross link sites identified by truncated cDNAs.

        [For reads with deletions] The deleted nucleotide in CLIP and iCLIP
        cDNAs mapped by Novoalign was used to define the cross-link sites
        identified by read-through cDNAs. If a cDNA had more than one deletion,
        we selected the one closest to the beginning of the read.

    The dtype to use internally for storage can be specified. Large types
    reduce the chance of overflow, but require more memory. With 'uint16'
    the largest count that can be handled is 255. Data is stored sparse,
    so memory is less of a problem. Overflow will cause a ValueError.

    returns a tuple of pandas Series objects, with the positive and negative
    strand arrays and also a counter object that contains the counts for each
    type of site."""

    pos_depths = collections.defaultdict(int)
    neg_depths = collections.defaultdict(int)

    counter = 0

    for read in reads:
        pos = getCrosslink(read)
        counter += 1

        if read.is_reverse:
            neg_depths[float(pos)] += 1
        else:
            pos_depths[float(pos)] += 1

    try:
        pos_depths = pd.Series(pos_depths, dtype=dtype)
    except ValueError:
        pos_depths = pd.Series({}, dtype=dtype)

    try:
        neg_depths = pd.Series(neg_depths, dtype=dtype)
    except ValueError:
        neg_depths = pd.Series({}, dtype=dtype)

    # check for integer overflow: counter sum should add up to array sum
    array_sum = pos_depths.sum() + neg_depths.sum()
    if not counter == array_sum:
        raise (
            ValueError,
            "Sum of depths is not equal to number of "
            "reads counted, possibly dtype %s not large enough" % dtype,
        )

    #    E.debug("Counted %i truncated on positive strand, %i on negative"
    #            % (counter.truncated_pos, counter.truncated_neg))
    #    E.debug("and %i deletion reads on positive strand, %i on negative"
    #            % (counter.deletion_pos, counter.deletion_neg))

    return (pos_depths, neg_depths, counter)
