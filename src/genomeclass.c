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

#include "genomeclass.h"
#include "alphabet.h"
#include "buffer.h"
#include "context.h"
#include "math.h"


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
long long int number_pos_data_sequence = 100000000;
int calculate_compression = 0;
int compression_geco = 0;
int compression_jarvis3 = 0;
int max_number_bases = 0;
int help_menu = 0;
int verbose = 0;
int calculate_entropy = 0;
int calculate_melting = 0;

pthread_mutex_t input_file_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t output_file_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t task_mutex = PTHREAD_MUTEX_INITIALIZER;

FILE *file;

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
    {"threads", required_argument, 0, 't'},
    {"verbose", no_argument, 0, 'v'},
    {0, 0, 0, 0}
};

// Print help menu
void program_usage(char *prog_path) {
    printf("\nUSAGE: .%s -t <number_of_threads> -i <input_fasta> -s -g -d <sequence_1> [sequence_n]...\n\n", strrchr(prog_path, '/'));
    printf("Program options -------------------------------------------------------------------------------\n");
    printf("-h, --help\t\tPrints this message\n");
    printf("-i, --input\t\tSet input file (FASTA format).\n");
    printf("-o, --output\t\tSet the output file (tsv format).\n");
    printf("-s, --size\t\tCalculates the size and the normalized size of the sequences.\n");
    printf("-g, --gc_content\tCalculates the GC content.\n");
    printf("-b, --base_percentage\tCalculates the percentage of the bases A, C, T, G and other in the sequence.\n");
    printf("-c, --compression\tCalculates the compressibility of the sequences (Markov models).\n");
    printf("-e, --entropy\t\tCalculates the entropy of the sequences.\n");
    printf("-m, --melting\t\tCalculates the maximum melting temperature.\n");
    printf("-x, --experiment\tCalculates the compressibility of the sequences (GeCo3).\n");
    printf("-j, --jarvis\t\tCalculates the compressibility of the sequences (JARVIS3).\n");
    printf("-d, --distance\t\tSet a sequence to calculate the distance (several sequences can be set).\n");
    printf("-t, --threads\t\tSets the number of threads.\n");
    printf("-v, --verbose\t\tVerbose mode - disables progress bar and prints the results.\n");

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
    while ((opt = getopt_long(argc, argv, "hi:o:sgbcemxjd:t:v", long_options, &option_index))  != -1) {
        
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

    int res = check_if_fa_or_fq(path_input_file, number_of_threads);

    if (res != 0) {
        exit(1);
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
 
    return 0;
}

// Initial passage over the input file, retrieves information regarding the sequences (single thread)
int initial_reading() {

    printf("\nStarting the analysis of the input file.\n");

    int ch;

    unsigned long long int index_file = 0;
    unsigned long long int index_data_sequences = -1;

    int is_header = 0;

    unsigned long long int length_seq = 0;
    unsigned long long int number_bases_seq = 0;

    unsigned long long int number_cg = 0;

    // Read each character until the end of file (EOF)
    while ((ch = fgetc(file)) != EOF) {

        if ((char)ch == '>') { // If the character is '>', then it is the begining of a new sequence


            if (index_data_sequences == -1) { // If it is the first sequence in a file, set index, else write the last position of the previous sequence

                index_data_sequences = 0;
                
            } else {

                data_all_sequences[index_data_sequences].end_sequence = index_file - 1;
                data_all_sequences[index_data_sequences].length_sequence = length_seq;
                data_all_sequences[index_data_sequences].number_bases = number_bases_seq;
                data_all_sequences[index_data_sequences].cg_content = number_cg;

                if (number_bases_seq > max_number_bases){ // Update longest sequence
                    max_number_bases = number_bases_seq;
                }

                index_data_sequences ++;
                length_seq = 0;

                // Reset count of the characters in the sequence
                data_all_sequences[index_data_sequences].number_a = 0;
                data_all_sequences[index_data_sequences].number_c = 0;
                data_all_sequences[index_data_sequences].number_t = 0;
                data_all_sequences[index_data_sequences].number_g = 0;
                data_all_sequences[index_data_sequences].number_other = 0;
                
                
            }


            if (index_data_sequences >= number_pos_data_sequence){ // If there are more than the pre defined number of sequences, increase the size of the array by 500 positions
    
                number_pos_data_sequence = number_pos_data_sequence + 500; // Increase the maximum number

                Seq_data *aux_array = malloc(number_pos_data_sequence * sizeof(Seq_data)); // Allocate the memory accordingly
                memset(aux_array, 0, number_pos_data_sequence * sizeof(Seq_data)); // Set memory to zeros

                for (int i = 0; i < index_data_sequences + 1; i++) { // Copy data to new array
                    aux_array[i] = data_all_sequences[i];
                }
                
                free(data_all_sequences); // Free the old array
                data_all_sequences = aux_array; // Assign the new array to the old name
                                
            }
            
            // Update the initial positions and set header
            data_all_sequences[index_data_sequences].init_header = index_file;
            is_header = 1;
        
        } else if (is_header == 1 && (char)ch == '\n') { // Middle/end of the header condition

            is_header = 0;
            data_all_sequences[index_data_sequences].end_header = index_file;
            length_seq = 0;
            number_bases_seq = 0;
            number_cg = 0;

        } else { // Is part of the sequence, increment sequence length and number of bases
            if ((char)ch != '\n') {
                number_bases_seq ++;

                if ((char)ch == 'c' || (char)ch == 'C' || (char)ch == 'g' || (char)ch == 'G' ) {
                    number_cg ++;
                }

                if ((char)ch == 'a' || (char)ch == 'A'){
                    data_all_sequences[index_data_sequences].number_a ++;

                } else if ((char)ch == 'c' || (char)ch == 'C'){
                    data_all_sequences[index_data_sequences].number_c ++;

                } else if ((char)ch == 't' || (char)ch == 'T'){
                    data_all_sequences[index_data_sequences].number_t ++;

                } else if ((char)ch == 'g' || (char)ch == 'G'){
                    data_all_sequences[index_data_sequences].number_g ++;

                } else {
                    data_all_sequences[index_data_sequences].number_other ++;
                }


            }

            length_seq ++;
        }

        index_file ++;
    }

    // Update the information of the last sequence
    data_all_sequences[index_data_sequences].end_sequence = index_file-1; 
    data_all_sequences[index_data_sequences].length_sequence = length_seq;
    data_all_sequences[index_data_sequences].number_bases = number_bases_seq;
    data_all_sequences[index_data_sequences].cg_content = number_cg;


    if (number_bases_seq > max_number_bases){
        max_number_bases = number_bases_seq;
    }

    printf("Max number bases %d\n\n", max_number_bases);

    number_sequences = index_data_sequences + 1;
    printf("Number of sequences in the file - %d\n", number_sequences);   
    
    return 0;

}

// Gets parts of a file given the start and end positions#include <stdio.h>
#include <stdlib.h>

int read_file_partially(unsigned long long int start_pos, unsigned long long int end_pos, char **content) {
    if (start_pos < 0 || end_pos < start_pos) {
        fprintf(stderr, "Invalid positions %d, %d\n", start_pos, end_pos);
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

    return 0;
}


// Calculates the distances between sequences set by the user (returns the probability of a given subsequence being the subseuqence set by the user)
Dist_Prob_sequence get_sequence_distance(char *content_sequence, char *subsequence, int number_bases_content_sequence){

    int number_times_subsequence_found = 0;
    int sum_distances = 0;
    int last_pos = 0;
    Dist_Prob_sequence results;
    int number_possibilities = number_bases_content_sequence - strlen(subsequence) + 1; 

    if (strlen(subsequence) < strlen(content_sequence)){ // If the subsequence is valid lengthwise
        
        for (int i = 0; i < strlen(content_sequence) - strlen(subsequence) + 1; i++) {

            if (content_sequence[i] == subsequence[0]){ // If current position is the begining of the subsequence, seek the rest

                int is_match = 1;

                for (int j = 1; j < strlen(subsequence); j++) { // Seeking the rest of the sequence

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


// Write the results to the output file (.tsv format)
int write_to_file(char* results){

    FILE *file_output;
    char column_name[150]; 

    if (access(output_path, F_OK) == 0) { // If the output file exists, append the results
        file_output = fopen(output_path, "a");
        if (!file_output) {
            perror("Failed to open the output file\n");
            return 1;
        }
    } else { // else, create the file and write the header
        file_output = fopen(output_path, "w"); 
        if (!file_output) {
            perror("Failed to open file\n");
            return 1;
        }

        char *first_line = malloc(sizeof(sequences_calc_distance) + 1000); // TODO - this may give some problems in the future depending on the number of sequences to seek distance   
        
        sprintf(first_line, "%s", "Sequence_id");  
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

        fprintf(file_output, "%s\n", first_line);  // Write the first line to the file
           
    }

    fprintf(file_output, "%s\n", results);  // Write the results to the file
    fclose(file_output); // Close the file

    return 0;

}

// Open file and get size
long long int get_size_file(char* file_name) {

    FILE *fp = fopen(file_name, "rb"); // open in binary mode
    fseek(fp, 0, SEEK_END);
    long size = ftell(fp);
    fclose(fp);

    return size;

}

// Create a file with a single sequence
void create_file_single_sequence(char * filename, char * sequence_to_write){

	if (access(filename, F_OK) != 0) { // If the file does not exist, create it
		// Open file for writing
		FILE * file_seq;
		file_seq = fopen(filename, "w");
		if (!file_seq) {
			perror("Failed to open uncompressed file");
			return -1;
		}

		// Write content
		fprintf(file_seq, "%s\n", sequence_to_write);
		fclose(file_seq);  // Flush and close before checking size
	}
}

float calculate_compression_ratio_geco (char * sequence_read, int id) {

    char filename_uncompressed[100];
    char filename_compressed[100];
    char logs_file[100];
    char command_geco[512];

    // Prepare filenames
    sprintf(filename_uncompressed, "sequence_%d.seq", id);
    sprintf(filename_compressed, "sequence_%d.seq.co", id);
    sprintf(logs_file, "logs_%d.txt", id);

	// Create the uncompressed file
	create_file_single_sequence(filename_uncompressed, sequence_read);

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

float calculate_compression_ratio_jarvis (char * sequence_read, int id) {

    char filename_uncompressed[100];
    char filename_compressed[100];
    char temp_fasta[100];
    char logs_file[100];
    char command_jarvis[512];
    char command_gto[512];
    char temp_name[512];

    // Prepare filenames
    sprintf(filename_uncompressed, "sequence_%d.seq", id);
    sprintf(temp_fasta, "sequence_%d.fa", id);
    sprintf(filename_compressed, "sequence_%d.seq.jc", id);
    sprintf(logs_file, "logs_%d_jarvis.txt", id);
    sprintf(temp_name, "tmp_%d.fa", id);

	// Create the uncompressed file
	FILE *fp = fopen(temp_fasta, "w");
    fprintf(fp, ">a\n%s\n", sequence_read);
    fclose(fp);

    // Build GTO command
    snprintf(command_gto, sizeof(command_gto), 
        "conda run -n genomeclass bash -c 'gto_fasta_rand_extra_chars < %s > tmp_%d.fa' && tail -n +2 tmp_%d.fa | tr -d '\n' > %s",
        temp_fasta, id, id, filename_uncompressed);

    

    // Run the command
    int ret = system(command_gto);
    if (ret != 0) {
        fprintf(stderr, "GTO command failed with return code %d\n", ret);
        return -1;
    }
    




    // Get size of uncompressed file
    long long int initial_size = get_size_file(filename_uncompressed);

    // Build JARVIS3 command
    snprintf(command_jarvis, sizeof(command_jarvis),
        "conda run -n genomeclass JARVIS3 -v -l 7 %s > %s 2>&1",
        filename_uncompressed, logs_file);

    // Run the command
    ret = system(command_jarvis);
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
    remove(temp_fasta);
    remove(filename_compressed);
    remove(temp_name);

    // Return compression ratio

    return (float) compressed_size / initial_size;

}

double calculate_compression_value(char * sequence_read, int id) {

    char name[100];
    int32_t ctx = 3;
	uint32_t alphaDen = 1;
	int32_t window_size = 2;
    int sym;
    double bits = 0;
    double ic = 0;
    uint64_t sequence_size = 0;
    uint8_t buf[BUFFER_SIZE];
    size_t bytes_read;

	// Prepare filenames
    sprintf(name, "sequence_%d.seq", id);

	create_file_single_sequence(name, sequence_read);

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

double calculate_entropy_value(char * sequence_read, int id) {

    char name[100];

    // Prepare filenames
    sprintf(name, "sequence_%d.seq", id);

    create_file_single_sequence(name, sequence_read);

    FILE *IN = fopen(name, "r");

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



// Tasks to be done by each thread
int worker_task(int index_data_sequence){

    char *content_header;
    char *content_sequence;

    //Initial and end position of the header
    unsigned long long int start_pos_header = data_all_sequences[index_data_sequence].init_header;
    unsigned long long int end_pos_header = data_all_sequences[index_data_sequence].end_header;

    // Get sequence header
    pthread_mutex_lock(&input_file_mutex);
    read_file_partially(start_pos_header, end_pos_header, &content_header);
    char* aux_header = content_header;
    pthread_mutex_unlock(&input_file_mutex);

    // Remove the \n and \t from the headers
    char *read_header = remove_newline_and_tab_characters(aux_header);

    //Initial and end position of the sequence   
    unsigned long long int start_pos_sequence = end_pos_header + 1;  // Start position in the file
    unsigned long long int end_pos_sequence = data_all_sequences[index_data_sequence].end_sequence;   // End position in the file

    // Get sequence
    pthread_mutex_lock(&input_file_mutex);
    read_file_partially(start_pos_sequence, end_pos_sequence, &content_sequence);
    char* aux_sequence = content_sequence;
    pthread_mutex_unlock(&input_file_mutex);

    // Remove \n characters from sequence
    char *read_sequence = remove_newline_and_tab_characters(aux_sequence);
    

    // Allocate the space required to store the results
    char *results = malloc(strlen(read_header) + sizeof(int) + sizeof(float) * (2 + number_sequences_calc_distance) + 100); // 100 more positions are assigned to avoid errors 
    
    // Copy the header to results
    sprintf(results, "%s",read_header);

    // Calculate the metrics
    if (calculate_size == 1) {
        int size_sequence = data_all_sequences[index_data_sequence].number_bases;
        float size_sequence_normalized = (float) data_all_sequences[index_data_sequence].number_bases / max_number_bases;

        // Copy results
        results = concatenate_strings(results, int_to_string(size_sequence), 1);
        results = concatenate_strings(results, float_to_string(size_sequence_normalized), 1);
    }

    if (calculate_gc_content == 1){
        float cg_content = (float) data_all_sequences[index_data_sequence].cg_content / data_all_sequences[index_data_sequence].number_bases;
        results = concatenate_strings(results, float_to_string(cg_content), 1);
    }

    if (calculate_base_percentage == 1) {
        float perc_A = (float) data_all_sequences[index_data_sequence].number_a / data_all_sequences[index_data_sequence].number_bases;
        float perc_C = (float) data_all_sequences[index_data_sequence].number_c / data_all_sequences[index_data_sequence].number_bases;
        float perc_T = (float) data_all_sequences[index_data_sequence].number_t / data_all_sequences[index_data_sequence].number_bases;
        float perc_G = (float) data_all_sequences[index_data_sequence].number_g / data_all_sequences[index_data_sequence].number_bases;
        float perc_Other = (float) data_all_sequences[index_data_sequence].number_other / data_all_sequences[index_data_sequence].number_bases;

        results = concatenate_strings(results, float_to_string(perc_A), 1);
        results = concatenate_strings(results, float_to_string(perc_C), 1);
        results = concatenate_strings(results, float_to_string(perc_T), 1);
        results = concatenate_strings(results, float_to_string(perc_G), 1);
        results = concatenate_strings(results, float_to_string(perc_Other), 1);
    }

    if (sequences_calc_distance != NULL) {
        
        // Copy results for each sub sequence considered
        for (int i = 0; i < number_sequences_calc_distance; i++){

            Dist_Prob_sequence sequence_data = get_sequence_distance(read_sequence, sequences_calc_distance[i], data_all_sequences[index_data_sequence].number_bases);
            float avg_sequence_distance = sequence_data.avg_distance;
            float sequence_probability = sequence_data.prob_sequence;
            
            results = concatenate_strings(results, float_to_string(avg_sequence_distance), 1);
            results = concatenate_strings(results, float_to_string(sequence_probability), 1);
        }
        
    }

    if (calculate_compression == 1) {
        double nc_results = calculate_compression_value(read_sequence, index_data_sequence);
        results = concatenate_strings(results, float_to_string(nc_results), 1);
    }

    if (calculate_entropy == 1) {
        double entropy_val = calculate_entropy_value(read_sequence, index_data_sequence);
        results = concatenate_strings(results, float_to_string(entropy_val), 1);
    }

    if (calculate_melting == 1) {
        double melting_temp = calculate_melting_temperature(data_all_sequences[index_data_sequence].number_a, data_all_sequences[index_data_sequence].number_t, data_all_sequences[index_data_sequence].number_c, data_all_sequences[index_data_sequence].number_g);
        results = concatenate_strings(results, float_to_string(melting_temp), 1);
    }
    
    if (compression_geco == 1) {
        float compression_ratio_geco = calculate_compression_ratio_geco(read_sequence, index_data_sequence);
        results = concatenate_strings(results, float_to_string(compression_ratio_geco), 1);

    }

    if (compression_jarvis3 == 1) {
        float compression_ratio_jarvis = calculate_compression_ratio_jarvis(read_sequence, index_data_sequence);
        results = concatenate_strings(results, float_to_string(compression_ratio_jarvis), 1);

    }

	char filename_to_delete[100];

    // Prepare filename
    sprintf(filename_to_delete, "sequence_%d.seq", index_data_sequence);
	if (access(filename_to_delete, F_OK) == 0) {
		remove(filename_to_delete);
	}

    // Write results to file
    pthread_mutex_lock(&output_file_mutex);
    write_to_file(results);
    pthread_mutex_unlock(&output_file_mutex);
    
    
    if (verbose == 0){ // Update the progress bar
        progress_bar(number_sequences);
    } else { // Print the results
        printf("%s\n", results);
    }
    return 0;

}

// Returns the index of the first sequence that isn't being worked on
int get_index_to_work_on () {

    if (number_tasks_assigned < number_sequences){
        pthread_mutex_lock(&task_mutex);
        int index = number_tasks_assigned;
        number_tasks_assigned ++;
        pthread_mutex_unlock(&task_mutex);
        return index;
    } else {
        return -1;
    }

}

// Thread function
void* thread_function(void* arg) {
    int index = get_index_to_work_on();

    while (index != -1) { // While there are tasks to be done, perform a task and ask for a new one

        worker_task(index);
        index = get_index_to_work_on();

    }

    return NULL;
}


int main(int argc, char *argv[]) {

    data_all_sequences = calloc(number_pos_data_sequence, sizeof(Seq_data));

    // Parse options and returns error if there's an issue
    int return_code = option_parsing(argc, argv);

    // If there is an error in the option parsing, exit
    if (return_code == 1) {
        printf("Error - Execution unsuccessful.\n");
        exit(1);  // Exit with error status
    } else if (help_menu == 1) {
        exit(0);
    }

    // If the input file does not exist, exit
    if (access(path_input_file, F_OK) != 0) {
        printf("Error - The input file does not exist.\n");
        exit(1);  // Exit with error status
    }

	// Open input file
    file = fopen(path_input_file, "r");
    if (file == NULL) {
        perror("Error opening file");
        return 1;
    }

    // Read the input file and get relevant info on the sequences
    initial_reading();

    pthread_t threads[number_of_threads];  // Array to hold thread IDs
    int thread_ids[number_of_threads];     // Array to hold thread arguments    
    
    // Initialize threads
    for (int i = 0; i < number_of_threads; i++) {
        thread_ids[i] = i;  // Assign unique ID to each thread
        if (pthread_create(&threads[i], NULL, thread_function, (void*)&thread_ids[i]) != 0) { //assign a function to the threads
            perror("Failed to create thread");
            exit(1);
        }
    }
    
    // Wait for all threads to finish
    for (int i = 0; i < number_of_threads; i++) {
        pthread_join(threads[i], NULL);
    }

    printf("\nAll threads have finished.\n");
	fclose(file);

}
