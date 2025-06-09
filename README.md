# GenomeClass
A tool for the analysis and classification of genome sequences in the FASTA format.

### REPLICATION ###

Run the whole experiment in a Linux system with:
<pre>
git clone https://github.com/mirakaya/GenomeClass
cd GenomeClass/src/
make clean
make
./genomeclass -i Samples/rand_500_sequences.fasta -s -g -c -e
</pre>

To see the possible execution options type
<pre>
./genomeclass -h
</pre>

This will print the following options:

```
USAGE: ./genomeclass -t <number_of_threads> -i <input_fasta> -s -g -d <sequence_1> [sequence_n]...

Program options -------------------------------------------------- ------------------------
-h, --help              Prints this message
-i, --input             Set input file (FASTA format).
-o, --output            Set the output file (tsv format).
-s, --size              Calculates the size and the normalized size of the sequences.
-g, --gc_content        Calculates the GC content.
-c, --compression       Calculates the compressibility of the sequences (Markov models).
-e, --entropy           Calculates the entropy of the seuqences.
-x, --experiment        Calculates the compressibility of the sequences (GeCo).
-d, --distance          Set a sequence to calculate the distance.
-t, --threads           Sets the number of threads.
-v, --verbose           Verbose mode - disables progress bar and prints the results.
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
