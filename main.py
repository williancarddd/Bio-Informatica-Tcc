
"""
Author: William Cardoso Barbosa
Github: williancarddd
e-mail: williancard123@gmail.com
""" 

"""
dimmer count: AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT
header: [legth_sequence,  base_a, base_t, base_c, base_g,  AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT, most_frequence_base_pattern, polindromes_count, palindrome_threshold, type_familie ]

RF03116 aCoV-5UTR
RF03117 bCoV-5UTR
RF03118 gCoV-5UTR
RF03119 dCoV-5UTR
RF03120 Sarbecovirus-5UTR (includes SARS-CoV-2)
RF03121 aCoV-3UTR
RF03122 bCoV-3UTR
RF03123 gCoV-3UTR
RF03124 dCoV-3UTR
RF03125 Sarbecovirus-3UTR (includes SARS-CoV-2)

Four families have been updated:

RF00164 Coronavirus s2m RNA
RF00165 Coronavirus 3’-UTR pseudoknot
RF00182 Coronavirus packaging signal 
RF00507 Coronavirus frameshifting stimulation element
"""

import csv, os

class Genetic:
  """
  Receipt [ATGTGCCGAGGCCACGCGGAGTACGATCGAGGGTACAGCACAA, ATATGCCGAGGCCACGCGGAGTACGATCGAGGGTACAGCATAA, ...., .]
  and transform the data in [value, value, value, value, value, value ]
  """

  def __init__(self, data, codeFilWithFamilie):
    self.data = data
    self.header = ['length_sequence', 'base_a', 'base_t', 'base_c', 'base_g',
     'AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT', 
     'most_frequence_base_pattern', 'polindromes_count', 'palindrome_threshold', 'type_familie']
    self.result = []
    self.codeFilWithFamilie = codeFilWithFamilie
    self.tetra_nucleotides = self._generate_tetratanucleotides()
    self.result.append(self.header)
    self._process_data()

  def _process_data(self):
    for sequence in self.data:
      self.result.append(self._process_sequence(sequence))

  def _process_sequence(self, sequence):
    """
    Process the sequence and return a list with the values
    """
    result = []
    result.append(len(sequence))
    result.append(sequence.count('A'))
    result.append(sequence.count('T'))
    result.append(sequence.count('C'))
    result.append(sequence.count('G'))
    result.extend(self._dimmer_count(sequence))
    result.append(self._get_most_frequence_base_pattern(sequence))
    result.append(self._get_polindromes_count(sequence))
    result.append(self._get_palindrome_threshold(sequence))
    result.append(self._get_famile_type())
    return result

  def _dimmer_count(self, sequence):
    """
    Return a list with the count of dimmers
    """
    dimmers = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    return [sequence.count(dimmer) for dimmer in dimmers]


  def _get_most_frequence_base_pattern(self, sequence):

    tetratanucleotides = self._generate_tetratanucleotides()
    count = [sequence.count(tetratanucleotides[i]) for i in range(len(tetratanucleotides))]
    return tetratanucleotides[count.index(max(count))]
    

  def _get_polindromes_count(self, sequence):
    polindromes = []
    for i in range(len(sequence)):
      for j in range(i, len(sequence)):
        if sequence[i:j] == sequence[i:j][::-1]:
          polindromes.append(sequence[i:j])
    return len(polindromes)

  def _get_palindrome_threshold(self, sequence):
    polindromes = []
    for i in range(len(sequence)):
      for j in range(i, len(sequence)):
        if sequence[i:j] == sequence[i:j][::-1]:
          polindromes.append(sequence[i:j])
    return len(max(polindromes))

  def _generate_tetratanucleotides(self):
    """
    Generate all tetratanucleotides
    """
    bases = ['A', 'T', 'C', 'G']
    tetratanucleotides = []
    for base1 in bases:
      for base2 in bases:
        for base3 in bases:
          for base4 in bases:
            tetratanucleotides.append(base1 + base2 + base3 + base4)
    return tetratanucleotides

  def _get_famile_type(self):
   
    if(self.codeFilWithFamilie == 'RF03116'):
      return 'alpha'
    elif(self.codeFilWithFamilie == 'RF03117'):
      return 'beta'
    elif(self.codeFilWithFamilie == 'RF03118'):
      return 'gamma'
    elif(self.codeFilWithFamilie == 'RF03119'):
      return 'delta'
    elif(self.codeFilWithFamilie == 'RF03120'):
      return 'Sarbecovirus-5UTR'
    elif(self.codeFilWithFamilie == 'RF03121'):
      return 'alpha'
    elif(self.codeFilWithFamilie == 'RF03122'):
      return 'beta'
    elif(self.codeFilWithFamilie == 'RF03123'):
      return 'gamma'
    elif(self.codeFilWithFamilie == 'RF03124'):
      return 'delta'
    elif(self.codeFilWithFamilie == 'RF03125'):
      return 'Sarbecovirus-3UTR'
    elif(self.codeFilWithFamilie == 'RF00164'):
      return 'Coronavirus s2m RNA'.replace(' ', '_')
    elif(self.codeFilWithFamilie == 'RF00165'):
      return "Coronavirus 3'-UTR pseudoknot".replace(' ', '_')
    elif(self.codeFilWithFamilie == 'RF00182'):
      return 'Coronavirus packaging signal'.replace(' ', '_')
    elif(self.codeFilWithFamilie == 'RF00507'):
      return 'Coronavirus frameshifting stimulation element'.replace(' ', '_')
    else:
      return 'Other'

  
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
  """
  /arffs
  /csv_result
  /fa_files
  """

  """
  file = ReadFileFa('RF00164.fa')
  data = file.read()
  genetic = Genetic(data)
  result = genetic.get_result()
  csv_writer = CsvWriter('result.csv', result[0])
  csv_writer.write(result[1:])
  """

  ##many files too
  for file_name in os.listdir('fa_files'):
    if file_name.endswith('.fa'):
      file_object = ReadFileFa('fa_files/' + file_name)
      data = file_object.read()
      genetic = Genetic(data, file_name.split('.')[0])
      result = genetic.get_result()
      csv_writer = CsvWriter('csv_result/' + file_name.replace('.fa', '.csv'), result[0])
      csv_writer.write(result[1:])
      print('file ' + file_name + ' done')
