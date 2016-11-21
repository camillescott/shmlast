from doit.task import clean_targets
import pandas as pd
import screed

from .profile import profile_task
from .util import create_doit_task as doit_task
from .util import ShortenedPythonAction, title, which


dna_to_aa={'TTT':'F','TTC':'F', 'TTA':'L','TTG':'L',
                'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
                'TAT':'Y','TAC':'Y', 'TAA':'X','TAG':'X','TGA':'X',
                'TGT':'C','TGC':'C', 'TGG':'W',
                'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
                'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
                'CAT':'H','CAC':'H', 'CAA':'Q','CAG':'Q',
                'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
                'ATT':'I','ATC':'I','ATA':'I', 'ATG':'M',
                'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
                'AAT':'N','AAC':'N', 'AAA':'K','AAG':'K',
                'AGT':'S','AGC':'S', 'AGA':'R','AGG':'R',
                'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
                'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
                'GAT':'D','GAC':'D', 'GAA':'E','GAG':'E',
                'GGT':'G','GGC':'G','GGA':'G','GGG':'G'}


__complementTranslation = { "A": "T", "C": "G", "G": "C", "T": "A", "N": "N" }
def complement(s):
    ''' Get the base complement.

    Args:
        s (str): The sequence to complement.
    '''
    c = "".join(__complementTranslation[n] for n in s)
    return c


def reverse(s):
    '''Reverse the sequence.

    Args:
        s (str): Sequence to reverse.
    '''
    r = "".join(reversed(s))
    return r


def peptides(seq, start):
    '''Translate the nucleotide sequence.

    Args:
        seq (str): The nucleotide sequence.
        start (int): Translation start position.
    Yields:
        A sequence of amino acid identifiers for each codon.
    '''

    for i in range(start, len(seq), 3):
        yield dna_to_aa.get(seq[i:i+3], "X")


def translate(seq):
    '''6-frame translation of the given nucleotide sequence.

    Args:
        seq (str): The nucleotide sequence.
    Yields:
        str: The translation in each frame.
    '''

    for i in range(3):
        pep = peptides(seq, i)
        yield "".join(pep)

    revcomp = reverse(complement((seq)))
    for i in range(3):
        pep = peptides(revcomp, i)
        yield "".join(pep)


def translate_fastx(input_fn, output_fn):
    '''Translate a nucleotide FASTA file.

    Args:
        input_fn (str): The FASTA file to translate.
        output_fn (str): Filename to store the results.
    '''

    with open(output_fn, 'w') as fp:
        for record in screed.open(input_fn):
            for frame, t in enumerate(translate(record.sequence)):
                name = '{0}_{1}'.format(record.name, frame)
                fp.write('>{0}\n{1}\n'.format(name, t))


@doit_task
@profile_task
def rename_task(input_fn, output_fn, name_map_fn='name_map.csv', prefix='tr'):
    '''Rename the FASTA idenfiers to play nicely with various programs.

    Args:
        input_fn (str): The FASTA to rename.
        output_fn (str): The filename of the renamed version.
        name_map_fn (str): Where to store the mapping of old to new names.
        prefix (str): Prefix to use for each transcript.
    Returns:
        dict: A doit task dictionary.
    '''
    
    def rename_input():
        name_map = []
        with open(output_fn, 'w') as output_fp:
            for n, record in enumerate(screed.open(input_fn)):
                new_name = '{0}{1}'.format(prefix, n)
                output_fp.write('>{0}\n{1}\n'.format(new_name,
                                                     record.sequence))
                name_map.append((record.name, new_name))

        pd.DataFrame(name_map,
                     columns=['old_name', 'new_name']).to_csv(name_map_fn,
                                                              index=False)

    return {'name': 'rename:{0}'.format(input_fn),
            'title': title,
            'actions': [ShortenedPythonAction(rename_input)],
            'targets': [output_fn, name_map_fn],
            'file_dep': [input_fn],
            'clean': [clean_targets]}


@doit_task
@profile_task
def translate_task(input_fn, output_fn): 
    '''Translate a nucleotide FASTA in six frames.

    Args:
        input_fn (str): The nucleotide FASTA.
        output_fn (str): Destination translated FASTA.
    Returns:
        dict: A doit task dictionary.
    '''

    return {'name': 'translate:{0}'.format(input_fn),
            'title': title,
            'actions': [ShortenedPythonAction(translate_fastx, args=[input_fn, output_fn])],
            'targets': [output_fn],
            'file_dep': [input_fn],
            'clean': [clean_targets]}

