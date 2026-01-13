#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

#include <time.h>
#include <getopt.h>
#include <pthread.h>
#include <string.h>
#include <errno.h>
#include <sys/ioctl.h>
#include <unistd.h>
#include <ctype.h>

#include "auxgenomeclass.h"
#include "alphabet.h"
#include "buffer.h"
#include "context.h"
#include "math.h"
#include "stats.h"


// Global variables
int number_of_threads = 1;
char *path_input_file = NULL; 
char **sequences_calc_distance = NULL; 
int number_sequences_calc_distance = 0;
char *output_path = "output.tsv";
int calculate_size = 0;
int calculate_gc_content = 0;
int calculate_base_percentage = 0;
Seq_data *data_all_sequences = NULL;
int number_sequences = 0;
int number_tasks_assigned = 0;
long long int number_pos_data_sequence = 10;
int calculate_compression = 0;
int compression_geco = 0;
int compression_jarvis3 = 0;
int max_number_bases = 0;
int help_menu = 0;
int verbose = 0;
int calculate_entropy = 0;
int calculate_melting = 0;
int calculate_additional_metrics = 0;

pthread_mutex_t output_file_mutex = PTHREAD_MUTEX_INITIALIZER;

int option_index = 0;

static struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"input", required_argument, 0, 'i'},
    {"output", required_argument, 0, 'o'},
    {"size", no_argument, 0, 's'},
    {"gc_content", no_argument, 0, 'g'},
    {"base_percentage", no_argument, 0, 'b'},
    {"normalized_compression", no_argument, 0, 'c'},
    {"entropy", no_argument, 0, 'e'},
    {"melting", no_argument, 0, 'm'},
    {"experiment", no_argument, 0, 'x'},
    {"jarvis", no_argument, 0, 'j'},
    {"distance", required_argument, 0, 'd'},
    {"additional_metrics", no_argument, 0, 'a'},
    {"threads", required_argument, 0, 't'},
    {"verbose", no_argument, 0, 'v'},
    {0, 0, 0, 0}
};

// Print help menu
void program_usage(char *prog_path) {
    printf("\nUSAGE: .%s -t <number_of_threads> -i <input_fasta> -s -g -d <sequence_1> [sequence_n]...\n\n", strrchr(prog_path, '/'));
    printf("Program options --------------------------------------------------------------------------------------------\n");
    printf("-h, --help\t\t\tPrints this message\n");
    printf("-i, --input\t\t\tSet input file (FASTA format).\n");
    printf("-o, --output\t\t\tSet the output file (tsv format).\n");
    printf("-s, --size\t\t\tCalculates the size and the normalized size of the sequences.\n");
    printf("-g, --gc_content\t\tCalculates the GC content.\n");
    printf("-b, --base_percentage\t\tCalculates the percentage of the bases A, C, T, G and other in the sequence.\n");
    printf("-c, --compression\t\tCalculates the compressibility of the sequences (Markov models).\n");
    printf("-e, --entropy\t\t\tCalculates the entropy of the sequences.\n");
    printf("-m, --melting\t\t\tCalculates the maximum melting temperature.\n");
    printf("-x, --experiment\t\tCalculates the compressibility of the sequences (GeCo3).\n");
    printf("-j, --jarvis\t\t\tCalculates the compressibility of the sequences (JARVIS3).\n");
    printf("-d, --distance\t\t\tSet a sequence to calculate the distance (several sequences can be set).\n");
    printf("-a, --additional_metrics\tCalculates additional metrics.\n");
    printf("-t, --threads\t\t\tSets the number of threads.\n");
    printf("-v, --verbose\t\t\tVerbose mode - disables progress bar and prints the results.\n");

    help_menu = 1;
}

// Gets the options selected by the user
int option_parsing(int argc, char *argv[]) {

    char *prog_path = argv[0];
    int opt;
    
    // If there are no arguments, show menu
    if ( argc <= 1) {
        program_usage(prog_path);
        return 0;
    } 

    // Input options
    while ((opt = getopt_long(argc, argv, "hi:o:sgbcemxjd:at:v", long_options, &option_index))  != -1) {
        
        switch (opt) {
            case 'h': 
                program_usage(prog_path);
                return 0;
            case 'i':
                path_input_file = optarg;
                if (path_input_file == NULL) {
                    printf("\nError: Input file not specified.\n");
                    return 1;
                }
                break;
            case 'o':
                output_path = optarg;
                break;
            case 's':
                calculate_size = 1;
                break;
            case 'g':
                calculate_gc_content = 1;
                break;
            case 'b':
                calculate_base_percentage = 1;
                break;
            case 'c':
                calculate_compression = 1;
                break;
            case 'e':
                calculate_entropy = 1;
                break;
            case 'm':
                calculate_melting = 1;
                break;
            case 'x':
                compression_geco = 1;
                break;
            case 'j':
                compression_jarvis3 = 1;
                break;
            case 'd':
                sequences_calc_distance = append(sequences_calc_distance, number_sequences_calc_distance, optarg);
                number_sequences_calc_distance ++;
                break;
            case 'a':
                calculate_additional_metrics = 1;
                break;
            case 't':
                number_of_threads = atoi(optarg);
                if (number_of_threads < 1) {
                    printf("The argument of option -t must be a positive integer.\n");
                    program_usage(prog_path);
                    return 1;
                }
                break;
            case 'v':
                verbose = 1;
                break;
            case '?':
                program_usage(prog_path);
                return 1;
        }
    }

    // If the output file exists, delete it
    if (access(output_path, F_OK) == 0) {
        remove(output_path);
    }


    // Print execution options
    printf("Number of threads: %d\n", number_of_threads);
    printf("Input file: %s\n", path_input_file ? path_input_file : "None");
    printf("Output file: %s\n", output_path ? output_path : "None");    

    if (number_sequences_calc_distance != 0) {
        printf("Sequences to check distances: ");
        for (int i = 0; i < number_sequences_calc_distance; i++) {
            printf("%s  ", sequences_calc_distance[i]);
        }
        printf("\n");

    } else {
        printf("Sequences to check distances: None\n");
    }

    if (path_input_file == NULL) {
        printf("No input file was selected, exiting...\n");
        return 1;
    }

    int res = check_if_fa_or_fq(path_input_file, number_of_threads);

    if (res != 0) {
        exit(1);
    }

 
    return 0;
}

// Gets parts of a file given the start and end positions
int read_file_partially(unsigned long long int start_pos, unsigned long long int end_pos, char **content) {
    
    if (end_pos < start_pos) {
        fprintf(stderr, "Invalid positions %llu, %llu\n", start_pos, end_pos);
        return 1;
    }

    // Calculate bytes to read once, ensure non-negative size
    size_t bytes_to_read = (size_t)(end_pos - start_pos + 1);

    // Allocate buffer before seeking to handle early errors
    char *buffer = malloc(bytes_to_read + 1);
    if (!buffer) {
        perror("Memory allocation failed");
        return 1;
    }

    // Open input file
    FILE *file = fopen(path_input_file, "r");
    if (file == NULL) {
        perror("Error opening file");
        return 1;
    }

    // Use fseeko instead of fseek for large files (optional, platform dependent)
    if (fseek(file, start_pos, SEEK_SET) != 0) {
        perror("Failed to seek file");
        free(buffer);
        return 1;
    }

    // fread may return less than requested, read in a loop until done
    size_t total_read = 0;
    while (total_read < bytes_to_read) {
        size_t read_now = fread(buffer + total_read, 1, bytes_to_read - total_read, file);
        if (read_now == 0) {
            if (feof(file)) break; // EOF reached early
            perror("Error reading file");
            free(buffer);
            return 1;
        }
        total_read += read_now;
    }

    // Null terminate the buffer
    buffer[total_read] = '\0';

    // Optionally shrink buffer if fewer bytes read (could realloc or just keep)
    // *content = realloc(buffer, total_read + 1); // optional, only if memory is a concern
    *content = buffer;
    fclose(file);

    return 0;
}


// Calculates the distances between sequences set by the user (returns the probability of a given subsequence being the subseuqence set by the user)
Dist_Prob_sequence get_sequence_distance(char *content_sequence, char *subsequence, int number_bases_content_sequence){

    int number_times_subsequence_found = 0;
    int sum_distances = 0;
    Dist_Prob_sequence results;
    int number_possibilities = number_bases_content_sequence - strlen(subsequence) + 1; 

    if (strlen(subsequence) < strlen(content_sequence)){ // If the subsequence is valid lengthwise
        
        for (size_t i = 0; i < strlen(content_sequence) - strlen(subsequence) + 1; i++) {

            if (content_sequence[i] == subsequence[0]){ // If current position is the begining of the subsequence, seek the rest

                int is_match = 1;

                for (size_t j = 1; j < strlen(subsequence); j++) {

                    if (content_sequence[i+j] != subsequence[j]) {
                        is_match = 0;
                    }

                }

                if (is_match == 1) { // If at the end of the matching it is still true, then add a new match to the sum_distances and number_of_sequences_found

                    if (number_times_subsequence_found != 0) {
                        sum_distances = i;
                    }

                    number_times_subsequence_found ++;
                    
                }

            } 

        }


        if (number_times_subsequence_found == 0){ // Error - no sequence found
            results.avg_distance = number_possibilities;
            results.prob_sequence = 0;
        } else {
            results.avg_distance = (float) sum_distances / number_times_subsequence_found;          
            results.prob_sequence = (float) number_times_subsequence_found / number_possibilities;
        }

        //printf("%f    %f\n", results.avg_distance , results.prob_sequence);
    
    } else { // Subsequence not valid lengthwise

        results.avg_distance = strlen(content_sequence);
        results.prob_sequence = 0;

    }

    return results;
    
}


// Write the results to an already opened output file (.tsv format)
int write_to_file(FILE *file_output, char* results) {

    if (!file_output) return 1;

    pthread_mutex_lock(&output_file_mutex);
    fprintf(file_output, "%s\n", results);  // Write the results to the file
    fflush(file_output); // Ensure it's written immediately
    pthread_mutex_unlock(&output_file_mutex);

    return 0;
}

// Write header line to output file
int write_header(FILE *file_output) {

    if (!file_output) return 1;

    char column_name[150]; 
    char *first_line = strdup("Sequence_id");  

    if (calculate_size == 1) {
        first_line = concatenate_strings(first_line, "Sequence_size", 1);
        first_line = concatenate_strings(first_line, "Normalized_sequence_size", 1);
    }
    if (calculate_gc_content == 1) {
        first_line = concatenate_strings(first_line, "CG_content", 1);
    }
    if (calculate_base_percentage == 1) {
        first_line = concatenate_strings(first_line, "Percentage_A\tPercentage_C\tPercentage_T\tPercentage_G\tPercentage_Other", 1);
    }
    if (sequences_calc_distance != NULL) {
        for (int i = 0; i < number_sequences_calc_distance; i++){
            sprintf(column_name, "Avg_distance_%s", sequences_calc_distance[i]);
            first_line = concatenate_strings(first_line, column_name, 1);
            sprintf(column_name, "Prob_sequence_%s", sequences_calc_distance[i]);
            first_line = concatenate_strings(first_line, column_name, 1);
        }
    }
    if (calculate_compression == 1) {
        first_line = concatenate_strings(first_line, "Compression_ratio(Markov_models)", 1);
    }
    if (calculate_entropy == 1){
        first_line = concatenate_strings(first_line, "Shannon_entropy", 1);
    }
    if (calculate_melting == 1){
        first_line = concatenate_strings(first_line, "Maximum_melting_temperature", 1);
    }
    if (compression_geco == 1) {
        first_line = concatenate_strings(first_line, "Compression_ratio(GeCo3)", 1);
    }
    if (compression_jarvis3 == 1) {
        first_line = concatenate_strings(first_line, "Compression_ratio(JARVIS3)", 1);
    }
    if (calculate_additional_metrics == 1) {
        first_line = concatenate_strings(first_line, 
            "Distinct_symbols\tRelative_diversity\tMost_common_symbol_freq_rel\tMax_min_freq_difference_rel\t"
            "Longest_run_length\tAverage_run_length\tChange_rate\tEqual_adjacent_pairs_ratio\t"
            "Distinct_bigrams\tBigram_diversity\tMost_frequent_bigram_1\tFreq_bigram_1\t"
            "Most_frequent_bigram_2\tFreq_bigram_2\tMost_frequent_bigram_3\tFreq_bigram_3", 1);
    }

    fprintf(file_output, "%s\n", first_line);
    fflush(file_output);
    free(first_line);

    return 0;
}


// Open file and get size
long long int get_size_file(const char *file_name) {
    FILE *fp = fopen(file_name, "rb"); // open in binary mode
    if (fp == NULL) {
        fprintf(stderr, "Error: could not open file '%s'\n", file_name);
        return -1;  // indicate failure
    }

    // Seek to end and get position
    if (fseek(fp, 0, SEEK_END) != 0) {
        fprintf(stderr, "Error: fseek failed for file '%s'\n", file_name);
        fclose(fp);
        return -1;
    }

    long long int size = ftell(fp);
    if (size == -1L) {
        fprintf(stderr, "Error: ftell failed for file '%s'\n", file_name);
        fclose(fp);
        return -1;
    }

    fclose(fp);
    return size;
}



// Removes all non-ACTG characters by replacing them with random valid bases
char *remove_non_base_chars(const char *sequence_to_write) {
    if (sequence_to_write == NULL) return NULL;

    const char valid_bases[] = "ACTG";
    int len = (int)strlen(sequence_to_write);

    // Allocate space for new sequence (+1 for null terminator)
    char *cleaned_seq = malloc(len + 1);
    if (!cleaned_seq) {
        perror("malloc failed");
        exit(EXIT_FAILURE);
    }

    srand((unsigned int) time(NULL)); // seed random generator

    for (int i = 0; i < len; i++) {
        char c = sequence_to_write[i];
        if (c == 'A' || c == 'C' || c == 'T' || c == 'G') {
            cleaned_seq[i] = c;
        } else if (c == 'a' || c == 'c' || c == 't' || c == 'g') {
            cleaned_seq[i] = (char) toupper(c);
        } else {
            cleaned_seq[i] = valid_bases[rand() % 4];
        }
    }

    cleaned_seq[len] = '\0';
    return cleaned_seq;
}

// Create a file with a single sequence
void create_file_single_sequence_seq(char * filename, char * sequence_to_write){

    //printf("Creating %s\n", filename);

    remove(filename);

    FILE * file_seq;
    file_seq = fopen(filename, "w");
    if (!file_seq) {
        perror("Failed to open uncompressed file");
        return;
    }

    char * clean_seq = remove_non_base_chars(sequence_to_write);

    // Write content
    fprintf(file_seq, "%s", clean_seq);
    fclose(file_seq);  // Flush and close before checking size

}

float calculate_compression_ratio_geco (int id, char* filename_uncompressed) {

    char filename_compressed[100];
    char logs_file[100];
    char command_geco[512];

    // Prepare filenames
    sprintf(filename_compressed, "sequence_%d.seq.co", id);
    sprintf(logs_file, "logs_%d.txt", id);

    // Get size of uncompressed file
    long long int initial_size = get_size_file(filename_uncompressed);

    // Build GeCo3 command
    snprintf(command_geco, sizeof(command_geco),
        "conda run -n genomeclass GeCo3 -tm 1:1:0:1:0.9/0:0:0 -tm 7:10:0:1:0/0:0:0 -tm 16:100:1:10:0/3:10:0.9 -lr 0.03 -hs 64 %s > %s 2>&1",
        filename_uncompressed, logs_file);

    // Run the command
    int ret = system(command_geco);
    if (ret != 0) {
        fprintf(stderr, "GeCo3 command failed with return code %d\n", ret);
        return -1;
    }

    // Get size of uncompressed file
    long long int compressed_size = get_size_file(filename_compressed);

    // Clean up files
    remove(filename_compressed);
    remove(logs_file);

    // Return compression ratio

    return (float) compressed_size / initial_size;

}

float calculate_compression_ratio_jarvis (int id, char *filename_uncompressed) {

    char filename_compressed[100];
    char logs_file[100];
    char command_jarvis[512];
    char temp_name[512];

    // Prepare filenames
    sprintf(filename_compressed, "sequence_%d.seq.jc", id);
    sprintf(logs_file, "logs_%d_jarvis.txt", id);
    sprintf(temp_name, "tmp_%d.fa", id);

    // Get size of uncompressed file
    long long int initial_size = get_size_file(filename_uncompressed);

    // Build JARVIS3 command
    snprintf(command_jarvis, sizeof(command_jarvis),
        "conda run -n genomeclass JARVIS3 -v -l 7 %s > %s 2>&1",
        filename_uncompressed, logs_file);

    // Run the command
    int ret = system(command_jarvis);
    if (ret != 0) {
        fprintf(stderr, "JARVIS3 command failed with return code %d\n", ret);
        return -1;
    }

    //printf("%s\n", filename_compressed);

    // Get size of uncompressed file
    long long int compressed_size = get_size_file(filename_compressed);

    // Clean up files
    remove(filename_compressed);
    remove(logs_file);
    remove(temp_name);

    // Return compression ratio

    return (float) compressed_size / initial_size;

}

double calculate_compression_value(int id, char* name) {

    int32_t ctx = 3;
	uint32_t alphaDen = 1;
	int32_t window_size = 2;
    int sym;
    double bits = 0;
    double ic = 0;
    uint64_t sequence_size = 0;
    uint8_t buf[BUFFER_SIZE];
    size_t bytes_read;

	FILE *IN = Fopen(name, "rb");

    ALPHABET *AL = CreateAlphabet();
    LoadAlphabet(AL, IN);

    //fprintf(stderr, "Alphabet cardinality: %u\n", AL->cardinality);

    CModel *CM = CreateCModel(ctx, alphaDen, AL->cardinality);
    CBUF   *symBuf = CreateCBuffer(BUFFER_SIZE, BGUARD);
    PModel *PM = CreatePModel(AL->cardinality);

    while((bytes_read = fread(buf, 1, BUFFER_SIZE, IN)) > 0) 
    for(size_t i = 0 ; i < bytes_read ; ++i) 
        {
        symBuf->buf[symBuf->idx] = sym = AL->revMap[buf[i]];
        GetPModelIdx(&symBuf->buf[symBuf->idx-1], CM);
        ComputePModel(CM, PM, CM->pModelIdx, CM->alphaDen);
        ic = PModelSymbolNats(PM, sym) / M_LN2;
        bits += ic;
        UpdateCModelCounter(CM, sym, CM->pModelIdx);
        UpdateCBuffer(symBuf);
        ++sequence_size;
        }

    RemovePModel(PM);
    RemoveCBuffer(symBuf);

    //fprintf(stderr, "NC: %lf\n", bits / ((double) sequence_size * log2(AL->cardinality)));

    fclose(IN);

    return bits / ((double) sequence_size * log2(AL->cardinality));

}


double calculate_entropy_value(int id, char* name) {

    FILE *IN = fopen(name, "r");
    if (IN == NULL) {
        fprintf(stderr, "ERROR: cannot open file %s\n", name);
        return -1.0;  // or handle cleanly
    }

    uint64_t freq[BYTE_RANGE] = {0};
    uint64_t total_bytes = 0;
    uint8_t  buffer[READ_BUFFER_SIZE];

    size_t bytes_read;
    while((bytes_read = fread(buffer, 1, READ_BUFFER_SIZE, IN)) > 0) 
    {
    total_bytes += bytes_read;
    for(size_t i = 0 ; i < bytes_read ; i++) 
        freq[buffer[i]]++;
    }

    fclose(IN);

    if(total_bytes == 0) 
    {
        fprintf(stderr, "Empty file.\n");
        exit(1);
    }

    double entropy = 0.0;
    for(int i = 0 ; i < BYTE_RANGE ; ++i) 
    {
    if(freq[i] == 0) continue;
    double p = (double)freq[i] / total_bytes;
    entropy -= p * log2(p);
    }

    //fprintf(stdout, "Shannon entropy: %.6f bits/byte\n", entropy);

    return entropy;


}

double calculate_melting_temperature (int number_A, int number_T, int number_C, int number_G) {

    uint32_t len = number_A + number_T + number_G + number_C;

    if(len == 0) return 0.0; // Avoid divide by zero

    if(len < 14) {
        return (number_A + number_T) * 2 + (number_G + number_C) * 4;
    } else {
        return 64.9 + 41.0 * (number_G + number_C - 16.4) / len;
    } 
    
}


char* make_results (int count_sequences, unsigned long long int start_pos_sequence, unsigned long long int current_pos, int number_a, int number_c, int number_t, int number_g, int number_other, char *header) {


    // Allocate the space required to store the results
    char *results = strdup("");

    int number_bases = number_a + number_c + number_g + number_t + number_other;
    float perc_A = (float) number_a / number_bases;
    float perc_C = (float) number_c / number_bases;
    float perc_T = (float) number_t / number_bases;
    float perc_G = (float) number_g / number_bases;
    float perc_Other = (float) number_other / number_bases;

    results = concatenate_strings(results, header, 0);    
    
    if (calculate_size == 1) {
        results = concatenate_strings(results, int_to_string(number_bases), 1);
        float normalized_size = (float) number_bases / max_number_bases;

        results = concatenate_strings(results, float_to_string(normalized_size), 1);
    }

    if (calculate_gc_content == 1){
        float cg_content = (float) (number_c + number_g ) / number_bases;
        results = concatenate_strings(results, float_to_string(cg_content), 1);
    }

    if (calculate_base_percentage == 1) {
        results = concatenate_strings(results, float_to_string(perc_A), 1);
        results = concatenate_strings(results, float_to_string(perc_C), 1);
        results = concatenate_strings(results, float_to_string(perc_T), 1);
        results = concatenate_strings(results, float_to_string(perc_G), 1);
        results = concatenate_strings(results, float_to_string(perc_Other), 1);
    }


    char seq_filename[100];                
    char *content_sequence;

    read_file_partially(start_pos_sequence, current_pos - 1, &content_sequence);
    char* aux_sequence = content_sequence;

    // Remove \n characters from sequence
    char *read_sequence = remove_newline_and_tab_characters(aux_sequence);

    // Create the seq file
    if (calculate_compression == 1 || compression_geco == 1 || compression_jarvis3 == 1 || calculate_entropy == 1 || calculate_additional_metrics == 1) {
        sprintf(seq_filename, "sequence_%d.seq", count_sequences);
        create_file_single_sequence_seq(seq_filename, read_sequence);
    }

    if (sequences_calc_distance != NULL) {

        // Copy results for each sub sequence considered
        for (int i = 0; i < number_sequences_calc_distance; i++){

            Dist_Prob_sequence sequence_data = get_sequence_distance(read_sequence, sequences_calc_distance[i], number_bases);
            float avg_sequence_distance = sequence_data.avg_distance;
            float sequence_probability = sequence_data.prob_sequence;
            
            results = concatenate_strings(results, float_to_string(avg_sequence_distance), 1);
            results = concatenate_strings(results, float_to_string(sequence_probability), 1);
        }
        
    }

    if (calculate_compression == 1) {
        double nc_results = calculate_compression_value(count_sequences, seq_filename);
        results = concatenate_strings(results, float_to_string(nc_results), 1);
    }

    if (calculate_entropy == 1) {
        double entropy_val = calculate_entropy_value(count_sequences, seq_filename);
        results = concatenate_strings(results, float_to_string(entropy_val), 1);
    }

    if (calculate_melting == 1) {
        double melting_temp = calculate_melting_temperature(number_a, number_t, number_c, number_g);
        results = concatenate_strings(results, float_to_string(melting_temp), 1);
    }
    
    if (compression_geco == 1) {
        float compression_ratio_geco = calculate_compression_ratio_geco(count_sequences, seq_filename);
        results = concatenate_strings(results, float_to_string(compression_ratio_geco), 1);

    }

    if (compression_jarvis3 == 1) {
        float compression_ratio_jarvis = calculate_compression_ratio_jarvis(count_sequences, seq_filename);
        results = concatenate_strings(results, float_to_string(compression_ratio_jarvis), 1);

    }

    if (calculate_additional_metrics == 1) {
        char * stats_results = stats_mode(seq_filename, 8);
        results = concatenate_strings(results, stats_results, 0);
    }

    remove(seq_filename);

    return results;
}

//Process the FASTA file and write output to TSV file
int process_file (int thread_id, FILE *file_output) {

    int ch;
    long int count_sequences = 0;
    int processing = 0;
    int is_header = 0;

    int number_a = 0;
    int number_c = 0;
    int number_t = 0;
    int number_g = 0;
    int number_other = 0;

    int max_header = 500;
    char header[max_header];
    int pos_header = 0;

    unsigned long long int start_pos_sequence = 0;

    unsigned long long int current_pos = 0;

    // Open input file
    FILE *file = fopen(path_input_file, "r");
    if (file == NULL) {
        perror("Error opening file");
        return 1;
    }



    // Read each character until the end of file (EOF)
    while ((ch = fgetc(file)) != EOF) {

        if ((char)ch == '>') { // If the character is '>', then it is the begining of a new sequence
            
            if (processing == 1) { // Was processing info; dump onto file

                char* results = make_results(count_sequences, start_pos_sequence, current_pos, number_a, number_c, number_t, number_g, number_other, header);
                
                // Write results to file
                write_to_file(file_output, results);

                // Update progress bar
                if (verbose == 0){ // Update the progress bar
                    progress_bar(number_sequences);
                    
                    //printf("Thread %d is processing sequence %ld\n", thread_id, count_sequences);
                } else { // Print the results
                    printf("%s\n", results);
                }
                pos_header = 0;
            }
            
            
            count_sequences ++;
            is_header = 1;

            // Reset variables for new sequence
            number_a = 0;
            number_c = 0;
            number_t = 0;
            number_g = 0;
            number_other = 0;

            if (count_sequences % number_of_threads == thread_id - 1) { //Process the content
                processing = 1;
                pos_header = 0;
                memset(header, 0, max_header);
                header[0] = '>';
                pos_header ++;
                header[pos_header] = '\0';

            } else {
                processing = 0;
                //printf("\nRejected %d to process sequence %ld, %d, res %ld\n", thread_id, count_sequences, number_of_threads, count_sequences % number_of_threads);

            }

            

        } else if (processing == 1) { // Either middle of a header or a part of the sequence

            if (is_header == 1) { // Then its a part of the header, copy content
                if ((char)ch == '\n') { // End of the header
                    is_header = 0;
                    start_pos_sequence = current_pos + 1;
                } else { // Copy content

                    if (pos_header < max_header) {
                        header[pos_header] = (char)ch;  
                        pos_header++;      
                        header[pos_header] = '\0';                   
                    }
                    
                }

            } else if ((char)ch != '\n') { // Is part of a sequence and is not a \n

                // Count chars
                if ((char)ch == 'a' || (char)ch == 'A'){
                    number_a ++;

                } else if ((char)ch == 'c' || (char)ch == 'C'){
                    number_c ++;

                } else if ((char)ch == 't' || (char)ch == 'T'){
                    number_t ++;

                } else if ((char)ch == 'g' || (char)ch == 'G'){
                    number_g ++;

                } else {
                    number_other ++;
                }


            }


        }

        current_pos ++;


    }

    // Process the last sequence in the file
    if (processing == 1) { // Was processing info; dump onto file
        char* results = make_results(count_sequences, start_pos_sequence, current_pos, number_a, number_c, number_t, number_g, number_other, header);

        // Write results to file
        write_to_file(file_output, results);

        // Update progress bar
        if (verbose == 0){ // Update the progress bar
            progress_bar(number_sequences);
            
            //printf("Thread %d is processing sequence %ld\n", thread_id, count_sequences);
        } else { // Print the results
            printf("%s\n", results);
        }
    }


    //printf("\nThread %d has finished processing.\n", thread_id);

    fclose(file);

    return 0;
}

// Thread function
void* thread_function(void* arg) {
    thread_arg_t *targ = (thread_arg_t*)arg;
    process_file(targ->thread_id, targ->file_output);
    return NULL;
}



int main(int argc, char *argv[]) {

    // Parse options and returns error if there's an issue
    int return_code = option_parsing(argc, argv);

    // If there is an error in the option parsing, exit
    if (return_code == 1) {
        //printf("Error - Execution unsuccessful.\n");
        exit(1);  // Exit with error status
    } else if (help_menu == 1) {
        exit(0);
    }

    // If the input file does not exist, exit
    if (access(path_input_file, F_OK) != 0) {
        printf("Error - The input file does not exist.\n");
        exit(1);  // Exit with error status
    }

    Info_file info_file = calculate_number_sequences(path_input_file);

    number_sequences = info_file.number_seqs;
    max_number_bases = info_file.max_size_seq;

    pthread_t threads[number_of_threads];  // Array to hold thread IDs
    int thread_ids[number_of_threads];     // Array to hold thread arguments 
    
    thread_arg_t thread_args[number_of_threads];
    
    FILE *file_output = fopen(output_path, "w");
    if (!file_output) {
        perror("Failed to open output file");
        exit(1);
    }

    write_header(file_output);
    
    // Initialize threads
    for (int i = 0; i < number_of_threads; i++) {
        thread_args[i].thread_id = i;
        thread_args[i].file_output = file_output;
        if (pthread_create(&threads[i], NULL, thread_function, (void*)&thread_args[i]) != 0) {
            perror("Failed to create thread");
            exit(1);
        }
    }
    
    // Wait for all threads to finish
    for (int i = 0; i < number_of_threads; i++) {
        pthread_join(threads[i], NULL);
    }

    fclose(file_output);

    printf("\nAll threads have finished.\n");

}