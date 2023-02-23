#! /usr/bin/env python3

import sys

def translate_sequence(rna_sequence, genetic_code):
    """Translates a sequence of RNA into a sequence of amino acids.

    Translates `rna_sequence` into string of amino acids, according to the
    `genetic_code` given as a dict. Translation begins at the first position of
    the `rna_sequence` and continues until the first stop codon is encountered
    or the end of `rna_sequence` is reached.

    If `rna_sequence` is less than 3 bases long, or starts with a stop codon,
    an empty string is returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    str
        A string of the translated amino acids.
    """
    if len(rna_sequence) < 3:
            return ""
    amino_acid_sequence = ""
    for i in range(0, len(rna_sequence), 3):
        codon = rna_sequence[i:i+3]
        amino_acid = genetic_code.get(codon.upper(),"")
        if amino_acid != "*":
            amino_acid_sequence += amino_acid
        else:
            break

    return amino_acid_sequence

def get_all_translations(rna_sequence, genetic_code):
    """Get a list of all amino acid sequences encoded by an RNA sequence.

    All three reading frames of `rna_sequence` are scanned from 'left' to
    'right', and the generation of a sequence of amino acids is started
    whenever the start codon 'AUG' is found. The `rna_sequence` is assumed to
    be in the correct orientation (i.e., no reverse and/or complement of the
    sequence is explored).

    The function returns a list of all possible amino acid sequences that
    are encoded by `rna_sequence`.

    If no amino acids can be translated from `rna_sequence`, an empty list is
    returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    list
        A list of strings; each string is an sequence of amino acids encoded by
        `rna_sequence`.
    """
    # Define list to store possible amino acid sequences
    amino_acid_sequences = []
    
    # Iterate over all possible reading frames
    for frame in range(3):
        # Extract subsequence for current reading frame
        rna_sequence=rna_sequence.upper()
        sub_sequence = rna_sequence[frame:]
        
        # Check if subsequence is long enough to encode an amino acid
        if len(sub_sequence) < 3:
            continue
        
        # Iterate over all codons in subsequence
        i = 0
        amino_acid_sequence = ""
        while i < len(sub_sequence) - 2:
            codon = sub_sequence[i:i+3]
            amino_acid = genetic_code.get(codon, None)
            
            # If current codon is a start codon, start a new amino acid sequence
            if amino_acid == "M":
                amino_acid_sequence = amino_acid
            
            # If amino acid sequence is not empty and current codon is a stop codon, add sequence to list
            elif amino_acid == "*":
                if len(amino_acid_sequence) > 0:
                    amino_acid_sequences.append(amino_acid_sequence)
                    amino_acid_sequence = ""
            
            # If amino acid sequence is not empty and current codon is not a stop codon, add amino acid to sequence
            elif len(amino_acid_sequence) > 0:
                amino_acid_sequence += amino_acid
            
            i += 3
        
        # Check if last codon in subsequence is a stop codon, and add sequence to list if necessary
        if amino_acid_sequence and amino_acid_sequence[-1] != "*":
            amino_acid_sequences.append(amino_acid_sequence)
    
    return amino_acid_sequences

def get_reverse(sequence):
    """Reverse orientation of `sequence`.
    Returns a string with `sequence` in the reverse order.
    If `sequence` is empty, an empty string is returned.
    Examples
    --------
    >>> get_reverse('AUGC')
    'CGUA'
    """
    reverse = sequence[::-1]
    reverse = reverse.upper()
    return reverse
    # if sequence:
    #    return sequence[::-1]
    #    else:
    #    print("")


def get_complement(sequence):
    """Get the complement of a `sequence` of nucleotides.
    Returns a string with the complementary sequence of `sequence`.
    If `sequence` is empty, an empty string is returned.
    Examples
    --------
    >>> get_complement('AUGC')
    'UACG'
    REVIEW
    """
    if  sequence:
        complement = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'}
        rna_complement = "".join([complement[base] for base in sequence.upper()])
        return rna_complement
    else:
         return ""


def reverse_and_complement(sequence):
    """Get the reversed and complemented form of a `sequence` of nucleotides.
    Returns a string that is the reversed and complemented sequence
    of `sequence`.
    If `sequence` is empty, an empty string is returned.
    Examples
    --------
    >>> reverse_and_complement('AUGC')
    'GCAU'
    REVIEW
    """
    if sequence:
        complement = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'}
        return ''.join(complement[base] for base in sequence.upper())[::-1]
    else:
            return ""

def get_longest_peptide(rna_sequence, genetic_code):
    """Get the longest peptide encoded by an RNA sequence.

    Explore six reading frames of `rna_sequence` (the three reading frames of
    `rna_sequence`, and the three reading frames of the reverse and complement
    of `rna_sequence`) and return (as a string) the longest sequence of amino
    acids that it encodes, according to the `genetic_code`.

    If no amino acids can be translated from `rna_sequence` nor its reverse and
    complement, an empty string is returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    str
        A string of the longest sequence of amino acids encoded by
        `rna_sequence`.
    """
    
    # Translate each reading frame and append them all into a list of sequences
    #Provide rna_sequence to sequence for use in the reverse_and_complement function
    AT=[]
    longest_peptide=''
    #Obtain reverse compliment of sequence
    RC=reverse_and_complement(rna_sequence)
    #get all translations of the rna-sequence
    ATRS=get_all_translations(rna_sequence, genetic_code)
    #redefine rna_sequence with the reverse compliment
    rna_sequence=RC
    #get all translations fo the reverse compliment
    ATRC= get_all_translations(rna_sequence, genetic_code)
    #Build list of "All Translations (AT)"
    AT=ATRC + ATRS
#Search total list for longest member
    if AT:
        longest_peptide = max(AT, key=len)
    return longest_peptide
    


if __name__ == '__main__':
    genetic_code = {'GUC': 'V', 'ACC': 'T', 'GUA': 'V', 'GUG': 'V', 'ACU': 'T', 'AAC': 'N', 'CCU': 'P', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'CAU': 'H', 'AAU': 'N', 'AGU': 'S', 'GUU': 'V', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'CCA': 'P', 'ACA': 'T', 'CCC': 'P', 'UGU': 'C', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'UGC': 'C', 'CAG': 'Q', 'GAU': 'D', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'AGG': 'R', 'GGG': 'G', 'UCC': 'S', 'UCA': 'S', 'UAA': '*', 'GGA': 'G', 'UAC': 'Y', 'GAC': 'D', 'UAG': '*', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'GGC': 'G', 'AUG': 'M', 'CUG': 'L', 'GAG': 'E', 'CUC': 'L', 'AGA': 'R', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'GCU': 'A', 'GAA': 'E', 'AUU': 'I', 'UUG': 'L', 'UUA': 'L', 'UGA': '*', 'UUC': 'F'}
    rna_seq = ("AUG"
            "UAC"
            "UGG"
            "CAC"
            "GCU"
            "ACU"
            "GCU"
            "CCA"
            "UAU"
            "ACU"
            "CAC"
            "CAG"
            "AAU"
            "AUC"
            "AGU"
            "ACA"
            "GCG")
    longest_peptide = get_longest_peptide(rna_sequence = rna_seq,
            genetic_code = genetic_code)
    assert isinstance(longest_peptide, str), "Oops: the longest peptide is {0}, not a string".format(longest_peptide)
    message = "The longest peptide encoded by\n\t'{0}'\nis\n\t'{1}'\n".format(
            rna_seq,
            longest_peptide)
    sys.stdout.write(message)
    if longest_peptide == "MYWHATAPYTHQNISTA":
        sys.stdout.write("Indeed.\n")
