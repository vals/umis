import regex as re
import itertools
from collections import defaultdict

def exact_barcode_filter(chunk, bc1, bc2):
    parser_re = re.compile('(.*):CELL_(?P<CB>.*):UMI_(.*)\\n(.*)\\n\\+\\n(.*)\\n')
    kept = []
    for read in chunk:
        match = parser_re.search(read).groupdict()
        cb1 = match['CB']
        if bc2:
            cb1, cb2 = cb1.split("-")
        if cb1 not in bc1:
            continue
        if bc2 and cb2 not in bc2:
            continue
        kept.append(read)
    return kept

def correcting_barcode_filter(chunk, bc1hash, bc2hash):
    parser_re = re.compile('(.*):CELL_(?P<CB>.*):UMI_(.*)\\n(.*)\\n\\+\\n(.*)\\n')
    kept = []
    for read in chunk:
        match = parser_re.search(read).groupdict()
        cb1 = match['CB']
        if bc2hash:
            cb1, cb2 = cb1.split("-")
        bc1corrected = bc1hash[cb1]
        if not bc1corrected:
            continue
        if bc2hash:
            bc2corrected = bc2hash[cb2]
            if not bc2corrected:
                continue
            correctbc = bc1corrected + "-" + bc2corrected
        else:
            correctbc = bc1corrected
        if correctbc == match['CB']:
            kept.append(read)
        else:
            read = read.replace("CELL_" + match['CB'], "CELL_" + correctbc)
            kept.append(read)
    return kept

def exact_sample_filter(read, barcodes):
    parser_re = re.compile('(.*):CELL_(.*):UMI_(.*):SAMPLE_(?P<SB>.*)\\n(.*)\\n\\+\\n(.*)\\n')
    match = parser_re.search(read).groupdict()
    sample = match['SB']
    if sample not in barcodes:
        return None
    return read

def correcting_sample_filter(read, barcodehash):
    parser_re = re.compile('(.*):CELL_(.*):UMI_(.*):SAMPLE_(?P<SB>.*)\\n(.*)\\n\\+\\n(.*)\\n')
    match = parser_re.search(read).groupdict()
    sample = match['SB']
    barcodecorrected = barcodehash[sample]
    if not barcodecorrected:
        return None
    correctbc = barcodecorrected
    if correctbc == match['SB']:
        return(read)
    else:
        read = read.replace("SAMPLE_" + match['SB'], "SAMPLE_" + correctbc)
        return(read)

class MutationHash(object):

    def __init__(self, strings, nedit):
        self.hash = mutationhash(strings, nedit)

    def __getitem__(self, barcode):
        result = self.hash[barcode]
        if len(result) != 1:
            return None
        else:
            return list(result)[0]

def mutationhash(strings, nedit):
    """
    produce a hash with each key a nedit distance substitution for a set of
    strings. values of the hash is the set of strings the substitution could
    have come from
    """
    maxlen = max([len(string) for string in strings])
    indexes = generate_idx(maxlen, nedit)
    muthash = defaultdict(set)
    for string in strings:
        muthash[string].update([string])
        for x in substitution_set(string, indexes):
            muthash[x].update([string])
    return muthash

def substitution_set(string, indexes):
    """
    for a string, return a set of all possible substitutions
    """
    strlen = len(string)
    return {mutate_string(string, x) for x in indexes if valid_substitution(strlen, x)}

def valid_substitution(strlen, index):
    """
    skip performing substitutions that are outside the bounds of the string
    """
    values = index[0]
    return all([strlen > i for i in values])

def generate_idx(maxlen, nedit):
    """
    generate all possible nedit edits of a string. each item has the form
    ((index1, index2), 'A', 'G')  for nedit=2
    index1 will be replaced by 'A', index2 by 'G'

    this covers all edits < nedit as well since some of the specified
    substitutions will not change the base
    """
    ALPHABET = ["A", "C", "G", "T"]
    indexlists = []
    ALPHABETS = [ALPHABET for x in range(nedit)]
    return list(itertools.product(itertools.combinations(range(maxlen), nedit),
                                  *ALPHABETS))

def mutate_string(string, tomutate):
    strlist = list(string)
    for i, idx in enumerate(tomutate[0]):
        strlist[idx] = tomutate[i+1]
    return "".join(strlist)
