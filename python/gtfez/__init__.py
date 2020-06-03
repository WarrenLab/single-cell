from collections import OrderedDict
import re


class GTFRecord:
    """
    A single line of a GTF file.

    Attributes:
        seqname: name of sequence (str)
        source: as in GTF specification (str)
        feature: as in GTF specification (str)
        start: as in GTF specification (int)
        end: as in GTF specification (int)
        score: as in GTF specification (str)
        strand: as in GTF specification (str)
        frame: as in GTF specification (str)
        attributes: dict mapping attribute names to their values
            in the attributes field of the GTF record
    """
    def __init__(self, record_string):
        """
        Makes a new GTFRecord from a single line of a GTF file.

        Args:
            record_string (str): a single line of a GTF file
        """
        fields = record_string.split('\t')
        if len(fields) < 8:
            raise Exception
        self.seqname, self.source, self.feature = fields[0:3]
        self.start, self.end = map(int, fields[3:5])
        self.score, self.strand, self.frame = fields[5:8]
        if len(fields) >= 9:
            self.attributes = parse_attributes(fields[8])

    def __str__(self):
        attributes_string = attributes_dict_to_string(self.attributes)
        return '\t'.join([
            self.seqname,
            self.source,
            self.feature,
            str(self.start),
            str(self.end),
            self.score,
            self.strand,
            self.frame,
            attributes_string,
        ])


attributes_re = re.compile('^(\w+) "(.+)"$')
def parse_attributes(attributes_string):
    """
    Parses the contents of a GTF attributes field into a dict.

    Args:
        attributes_string: The attributes string of a GTF record;
            i.e., the 9th field

    Returns:
        an OrderedDict mapping attribute keywords to their values

    >>> parse_attributes('gene_name "ESR1"; gene_biotype "protein_coding";')
    OrderedDict([('gene_name', 'ESR1'), ('gene_biotype', 'protein_coding')])
    """
    # break attributes into individual 'key "val"' strings, removing
    # whitespace and the empty one after the last ';'
    attribute_strings = map(str.strip, attributes_string.split(';')[:-1])

    # read all of the pairs into a dictionary and return
    attributes_dict = OrderedDict()
    for s in attribute_strings:
        match = attributes_re.match(s)
        attributes_dict[match.group(1)] = match.group(2)
    return attributes_dict


def attributes_dict_to_string(attributes_dict):
    """
    Converts an attributes dict back into GTF string format.

    Args:
        attributes_dict: a dict mapping attribute keywords to their
            values

    Returns:
        a string of the given attributes in GTF format

    >>> attributes_dict_to_string({
    ...     'gene_name': 'ESR1',
    ...     'gene_biotype': 'protein_coding',
    ... })
    'gene_name "ESR1"; gene_biotype "protein_coding";'
    """
    output_strings = []
    for key, value in attributes_dict.items():
        output_strings.append('{} "{}";'.format(key, value))
    return ' '.join(output_strings)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
