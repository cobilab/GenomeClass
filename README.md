# Classification tool
A tool for the classification of genome sequences in the FASTA format.

### REPLICATION ###

Run the whole experiment in a Linux system with:
<pre>
git clone https://github.com/mirakaya/classification_tool_C
gcc -O3 ct.c -pthread
./a.out -i Samples/a.fasta -t 3
</pre>

To see the possible reconstruction options type
<pre>
./a.out -h
</pre>

This will print the following options:

```
USAGE: ./a.out -t <number_of_threads> -i <input_fasta> -s -g -d <sequence_1> [sequence_n]...

Program options -------------------------------------------------- ------------------------
-h, --help              Prints this message
-i, --input             Set input file (FASTA format).
-o, --output            Set the output file (tsv format).
-s, --size              Calculates the size and the normalized size of the sequences.
-g, --gc_content        Calculates the GC content.
-c, --compression       Calculates the compressibility of the sequences (Markov models).
-x, --experiment        Calculates the compressibility of the sequences (GeCo).
-d, --distance          Set a sequence to calculate the distance.
-t, --threads           Sets the number of threads.
-v, --verbose           Verbose mode - disables progress bar and prints the results.
```

### CITATION ###

On using this software/method please cite:

* pending

### ISSUES ###

For any issue let us know at [issues link](https://github.com/mirakaya/classification_tool_C/issues).

### LICENSE ###

GPL v3.

For more information:
<pre>http://www.gnu.org/licenses/gpl-3.0.html</pre>
