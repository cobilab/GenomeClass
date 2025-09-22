# GenomeClass
A tool for the analysis and classification of genome sequences in the FASTA format.

### REPLICATION ###

Install the necessary dependencies in a Linux system with:
<pre>
git clone https://github.com/mirakaya/GenomeClass
cd GenomeClass/src/
./Installation.sh
</pre>

To download the datasets type
<pre>
./Download_datasets.sh
</pre>

To analyse a dataset type
<pre>
make clean
make
./genomeclass -i sequences.fasta -s -g -c -e -m -t 4
</pre>

To see the possible analysis options type
<pre>
./genomeclass -h
</pre>

This will print the following options:

```
USAGE: ./genomeclass -t <number_of_threads> -i <input_fasta> -s -g -d <sequence_1> [sequence_n]...

Program options -------------------------------------------------------------------------------
-h, --help              Prints this message
-i, --input             Set input file (FASTA format).
-o, --output            Set the output file (tsv format).
-s, --size              Calculates the size and the normalized size of the sequences.
-g, --gc_content        Calculates the GC content.
-b, --base_percentage   Calculates the percentage of the bases A, C, T, G and other in the sequence.
-c, --compression       Calculates the compressibility of the sequences (Markov models).
-e, --entropy           Calculates the entropy of the sequences.
-m, --melting           Calculates the maximum melting temperature.
-x, --experiment        Calculates the compressibility of the sequences (GeCo3).
-j, --jarvis            Calculates the compressibility of the sequences (JARVIS3).
-d, --distance          Set a sequence to calculate the distance (several sequences can be set).
-t, --threads           Sets the number of threads.
-v, --verbose           Verbose mode - disables progress bar and prints the results.
```

To train the ML models and classify the sequences in a dataset type
<pre>
python3 genomeclass.py -tf sequences_train.fasta -cf sequences_to_classify.fasta
</pre>

To see the possible classification options type
<pre>
python3 genomeclass.py -h
</pre>

This will print the following options:
```
usage: genomeclass.py [-h] [-tf <training_file>] [-tt <training_file>]
                      [-cf <file_to_classify>] [-ct <file_to_classify>]
                      [-s <position>] [-o <analysis_options>]
                      [-p <number_bases_permutations>] [-b]

Index

options:
  -h, --help            show this help message and exit
  -tf <training_file>, --training_fasta <training_file>
                        Input training multi-FASTA file
  -tt <training_file>, --training_tsv <training_file>
                        Input training TSV file
  -cf <file_to_classify>, --classification_fasta <file_to_classify>
                        Input FASTA file containing the sequences to be
                        classified
  -ct <file_to_classify>, --classification_tsv <file_to_classify>
                        Input TSV file containing the sequences to be
                        classified
  -s <position>, --segment <position>
                        Part of the Sequence_id that will become the target
                        feature
  -o <analysis_options>, --analysis_options <analysis_options>
                        Options for the execution of the C file. Please
                        surround the options with ""
  -p <number_bases_permutations>, --permutations <number_bases_permutations>
                        Add permutations of a certain number of characters
  -b, --balance         Balances the training dataset

```


### CITATION ###

On using this software/method please cite:

* pending

### ISSUES ###

For any issue let us know at [issues link](https://github.com/mirakaya/GenomeClass/issues).

### LICENSE ###

GPL v3.

For more information:
<pre>http://www.gnu.org/licenses/gpl-3.0.html</pre>
