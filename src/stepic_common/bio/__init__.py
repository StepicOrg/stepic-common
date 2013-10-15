import random

import Bio.Data.IUPACData

from ..common import rand_N


def generate_DNA(N):
    return ''.join(random.choice('ACGT') for _ in range(N))


def generate_protein(N):
    return ''.join(random.choice(Bio.Data.IUPACData.protein_letters) for _ in range(N))


def generate_repetitive_DNA(N):
    parts = 5
    part_len = min(30, N // 100) + 1
    parts = [generate_DNA(rand_N(part_len)) for _ in range(parts)]
    dna = []
    dna_len = 0
    while dna_len + part_len < N:
        part = random.choice(parts)
        dna.append(part)
        dna_len += len(part)
    return ''.join(dna)


def rc(s):
    return s.translate('*****************************************************************TVGHEFCDIJMLKNOPQYSAUBWXRZ[\]^_`tvghefcdijmlknopqysaubwxrz*************************************************************************************************************************************')[::-1]
