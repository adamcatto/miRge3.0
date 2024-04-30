from cutadapt.modifiers import SingleEndModifier
from dnaio import SequenceRecord
import re


class RepeatNucleotideTrimmer(SingleEndModifier):
    """Trims Ns repeated ≥ n times from the 3' and 5' end of reads"""

    def __init__(self, n=7):
        pattern = r'(A{%d,}|C{%d,}|G{%d,}|T{%d,})' % (n, n, n, n)
        self.regex = re.compile(pattern)

    def __repr__(self):
        return "RepeatNucleotideTrimmer()"

    def __call__(self, read, info):
        # Find all matches of the repeat pattern in the sequence
        matches = [(match.start(), match.end()) for match in self.regex.finditer(str(read.sequence))]

        # Initialize lists to store the trimmed sequence and corresponding qualities
        trimmed_sequence = []
        trimmed_qualities = []

        # Iterate over each base in the original sequence and check if it is retained after trimming
        current_pos = 0
        for start, end in matches:
            # Add non-trimmed bases and their corresponding qualities to the trimmed lists
            trimmed_sequence.extend(read.sequence[current_pos:start])
            trimmed_qualities.extend(read.qualities[current_pos:start])
            # Update current position to the end of the trimmed region
            current_pos = end

        # Add remaining non-trimmed bases and their corresponding qualities to the trimmed lists
        trimmed_sequence.extend(read.sequence[current_pos:])
        trimmed_qualities.extend(read.qualities[current_pos:])

        # Create a new SequenceRecord with the trimmed sequence and adjusted qualities
        trimmed_read = SequenceRecord(
            read.name, 
            ''.join(trimmed_sequence), 
            ''.join(trimmed_qualities)
        )

        return trimmed_read