"""PhiX Indexed Control (1000 Cycle) from Illumina

This module contains the index sequences for PhiX Indexed Control (1000 Cycle),
a spike-in control designed for high-cycle sequencing and low-diversity library pools.

Source: https://knowledge.illumina.com/instrumentation/general/instrumentation-general-reference_material-list/000009867

Key Features:
- 5 unique dual-index combinations (10 bp each)
- Total library size: ~825 bp (700 bp insert)
- Designed for up to 1000 cycles
- Provides color balancing during index reads
- Index combinations do not overlap with current Illumina or IDT indexes

Usage:
Apply for "High-Risk" runs:
- Low index diversity pools (< 5 unique index combinations)
- Long-read chemistry (> 350 total cycles)
"""

# PhiX Indexed Control index combinations
# Format: (i7_sequence, i5_forward_sequence)
PhiX_indexed_control = {
    "PhiX-Index1": {
        "i7": "CGACCTAACG",
        "i5": "ATACGCCGGC",
        "name": "PhiX-Index1",
        "num": 1,
    },
    "PhiX-Index2": {
        "i7": "TAGTTCGGTA",
        "i5": "CTGTCACTTA",
        "name": "PhiX-Index2",
        "num": 2,
    },
    "PhiX-Index3": {
        "i7": "ACCGGCCGTA",
        "i5": "GGCGCCATTG",
        "name": "PhiX-Index3",
        "num": 3,
    },
    "PhiX-Index4": {
        "i7": "AATCACCAGC",
        "i5": "CAGTAATTAC",
        "name": "PhiX-Index4",
        "num": 4,
    },
    "PhiX-Index5": {
        "i7": "AGCGAATTAG",
        "i5": "GGCGATCAGA",
        "name": "PhiX-Index5",
        "num": 5,
    },
}

# List of all PhiX index pairs for easy iteration
PhiX_index_pairs = [
    ("CGACCTAACG", "ATACGCCGGC"),  # PhiX-Index1
    ("TAGTTCGGTA", "CTGTCACTTA"),  # PhiX-Index2
    ("ACCGGCCGTA", "GGCGCCATTG"),  # PhiX-Index3
    ("AATCACCAGC", "CAGTAATTAC"),  # PhiX-Index4
    ("AGCGAATTAG", "GGCGATCAGA"),  # PhiX-Index5
]

# Helper function to get all PhiX indexes as tuples
def get_phix_index_pairs():
    """Return all PhiX index pairs as list of (i7, i5) tuples."""
    return PhiX_index_pairs


# Helper function to check if an index pair matches PhiX
def is_phix_index(i7_seq, i5_seq):
    """Check if given index pair matches any PhiX index combination."""
    return (i7_seq.upper(), i5_seq.upper()) in [
        (pair[0].upper(), pair[1].upper()) for pair in PhiX_index_pairs
    ]
