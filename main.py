
"""
Author: William Cardoso Barbosa
Github: williancarddd
""" 

"""
header: [legnth_sequence, base_count, pair_bases, most_frequence_base_pattern, polindromes_count, palindrome_threshold ]
"""

import csv

class Genetic:
  """
  Receipt [ATGTGCCGAGGCCACGCGGAGTACGATCGAGGGTACAGCACAA, ATATGCCGAGGCCACGCGGAGTACGATCGAGGGTACAGCATAA, ...., .]
  and transform the data in [value, value, value, value, value, value ]
  """

  def __init__(self, data):
    self.data = data
    self.header = ['length_sequence', 'base_count', 'pair_bases', 'most_frequence_base_pattern', 'polindromes_count', 'palindrome_threshold']
    self.result = []
    self.result.append(self.header)
    self._process_data()

  def _process_data(self):
    for sequence in self.data:
      self.result.append(self._process_sequence(sequence))

  def _process_sequence(self, sequence):
    length_sequence = len(sequence)
    pair_bases = self._get_pair_bases(sequence)
    base_count = len(pair_bases)
    most_frequence_base_pattern = self._get_most_frequence_base_pattern(sequence)
    polindromes_count = self._get_polindromes_count(sequence)
    palindrome_threshold = self._get_palindrome_threshold(sequence)
    return [length_sequence, base_count, pair_bases, most_frequence_base_pattern, polindromes_count, palindrome_threshold]

  def _get_pair_bases(self, sequence):
    bases = ['A', 'T', 'C', 'G']
    pairs = []
    for base in bases:
      pairs.append(sequence.count(base))
    return pairs

  def _get_most_frequence_base_pattern(self, sequence):
    bases = ['A', 'T', 'C', 'G']
    pairs = self._get_pair_bases(sequence)
    return bases[pairs.index(max(pairs))]

  def _get_polindromes_count(self, sequence):
    polindromes = []
    for i in range(len(sequence)):
      for j in range(i, len(sequence)):
        if sequence[i:j] == sequence[i:j][::-1]:
          polindromes.append(sequence[i:j])
    return len(polindromes)

  def _get_palindrome_threshold(self, sequence):
    """
    Palindrome Threshold returns the length of the palindrome which is maximum.

    """
    polindromes = []
    for i in range(len(sequence)):
      for j in range(i, len(sequence)):
        if sequence[i:j] == sequence[i:j][::-1]:
          polindromes.append(sequence[i:j])
    return len(max(polindromes))

  def get_result(self):
    return self.result

class ReadFileFa:
  """
  recept a:
    >FJ376620.1/26255-26297 Bulbul coronavirus HKU11-796, complete genome.
    ATGTGCCGAGGCCACGCGGAGTACGATCGAGGGTACAGCACAA
    >FJ376621.1/26173-26215 Thrush coronavirus HKU12-600, complete genome.
    ATATGCCGAGGCCACGCGGAGTACGATCGAGGGTACAGCATAA
    >.....
    ....
    ...
    .
    generate a list of sequence with [ATGTGCCGAGGCCACGCGGAGTACGATCGAGGGTACAGCACAA, ATATGCCGAGGCCACGCGGAGTACGATCGAGGGTACAGCATAA, ...., .]
  """
  def __init__(self, file):
      self.file = file
      self.sequence = []

  def read(self):
    with open(self.file, 'r') as f:
      for line in f:
        if line[0] != '>':
          self.sequence.append(line.strip())

    return self.sequence

class CsvWriter:
    def __init__(self, file_name, header):
        self.file_name = file_name
        self.header = header

    def write(self, data):
        with open(self.file_name, 'w') as file:
            writer = csv.writer(file)
            writer.writerow(self.header)
            writer.writerows(data)

if __name__ == '__main__':
  file = ReadFileFa('RF00164.fa')
  data = file.read()
  genetic = Genetic(data)
  result = genetic.get_result()
  csv_writer = CsvWriter('result.csv', result[0])
  csv_writer.write(result[1:])

