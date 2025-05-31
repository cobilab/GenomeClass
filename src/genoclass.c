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

#include "genoclass.h"


// Global variables
int number_of_threads = 1;
char *path_input_file = NULL; 
char **sequences_calc_distance = NULL; 
int number_sequences_calc_distance = 0;
char *output_path = "output.tsv";
int calculate_size = 0;
int calculate_gc_content = 0;
Seq_data *data_all_sequences = NULL;
int number_sequences = 0;
int number_tasks_assigned = 0;
int number_pos_data_sequence = 10;
int calculate_compression = 0;
int compression_geco;
int max_number_bases = 0;
int help_menu = 0;
int verbose = 0;

pthread_mutex_t input_file_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t output_file_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t task_mutex = PTHREAD_MUTEX_INITIALIZER;

int option_index = 0;

static struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"input", required_argument, 0, 'i'},
    {"output", required_argument, 0, 'o'},
    {"size", no_argument, 0, 's'},
    {"gc_content", no_argument, 0, 'g'},
    {"compression", no_argument, 0, 'c'},
    {"experiment", no_argument, 0, 'x'},
    {"distance", required_argument, 0, 'd'},
    {"threads", required_argument, 0, 't'},
    {"verbose", no_argument, 0, 'v'}
};

// Print help menu
void program_usage(char *prog_path) {
    printf("\nUSAGE: .%s -t <number_of_threads> -i <input_fasta> -s -g -d <sequence_1> [sequence_n]...\n\n", strrchr(prog_path, '/'));
    printf("Program options -------------------------------------------------- ------------------------\n");
    printf("-h, --help\t\tPrints this message\n");
    printf("-i, --input\t\tSet input file (FASTA format).\n");
    printf("-o, --output\t\tSet the output file (tsv format).\n");
    printf("-s, --size\t\tCalculates the size and the normalized size of the sequences.\n");
    printf("-g, --gc_content\tCalculates the GC content.\n");
    printf("-c, --compression\tCalculates the compressibility of the sequences (Markov models).\n");
    printf("-x, --experiment\tCalculates the compressibility of the sequences (GeCo).\n");
    printf("-d, --distance\t\tSet a sequence to calculate the distance.\n");
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
    while ((opt = getopt_long(argc, argv, "hi:o:sgcxd:t:v", long_options, &option_index))  != -1) {
        
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
            case 'c':
                calculate_compression = 1;
                break;
            case 'x':
                compression_geco = 1;
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

    FILE *file;
    int ch;

    int index_file = 0;
    int index_data_sequences = -1;

    int is_header = 0;

    int length_seq = 0;
    int number_bases_seq = 0;

    int number_cg = 0;

    // Open input file
    file = fopen(path_input_file, "r");
    if (file == NULL) {
        perror("Error opening file");
        return 1;
    }

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

    fclose(file);  // Close the file

    number_sequences = index_data_sequences + 1;
    printf("Number of sequences in the file - %d\n", number_sequences);   
    
    return 0;

}

// Gets parts of a file given the start and end positions
int read_file_partially (int start_pos, int end_pos, char *file_name, char **content) {

    // Open file
    FILE *file = fopen(file_name, "r");
    if (file == NULL) {
        perror("Failed to open file");
        return 1;
    }

    // Ensure the start and end positions are valid
    if (start_pos < 0 || end_pos < start_pos) {
        fprintf(stderr, "Invalid positions %d, %d\n", start_pos, end_pos);
        fclose(file);
        return 1;
    }

    // Move to the start position
    fseek(file, start_pos, SEEK_SET);

    // Calculate how many bytes to read
    int bytes_to_read = end_pos - start_pos + 1;  // +1 to include the end position

    // Allocate buffer to store the content
    char *buffer = malloc(bytes_to_read + 1);  // +1 for null terminator
    if (buffer == NULL) {
        perror("Memory allocation failed");
        fclose(file);
        return 1;
    }

    // Read the content into the buffer
    size_t bytes_read = fread(buffer, 1, bytes_to_read, file);
    if (bytes_read != bytes_to_read) {
        perror("Error reading the file");
        free(buffer);
        fclose(file);
        return 1;
    }

    // Null-terminate the buffer to make it a string
    buffer[bytes_read] = '\0';

    // Output the buffer read
    *content = buffer;

    // Close file
    fclose(file);

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

    FILE *file;
    char column_name[150]; 

    if (access(output_path, F_OK) == 0) { // If the output file exists, append the results
        file = fopen(output_path, "a");
        if (!file) {
            perror("Failed to open file");
            return 1;
        }
    } else { // else, create the file and write the header
        file = fopen(output_path, "w"); 
        if (!file) {
            perror("Failed to open file");
            return 1;
        }

        char *first_line = malloc(sizeof(sequences_calc_distance) + 1000); // TODO - this may give some problems in the future depending on the number of sequences to seek distance   
        
        sprintf(first_line, "%s", "Sequence id");  
        if (calculate_size == 1) {
            first_line = concatenate_strings(first_line, "Sequence size", 1);
            first_line = concatenate_strings(first_line, "Normalized sequence size", 1);
        }
        if (calculate_gc_content == 1) {
            first_line = concatenate_strings(first_line, "CG content", 1);
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
            first_line = concatenate_strings(first_line, "Compression rate (Markov models)", 1);
        }

        if (compression_geco == 1) {
            first_line = concatenate_strings(first_line, "Compression rate (GeCo3)", 1);
        }

        fprintf(file, "%s\n", first_line);  // Write the first line to the file
           
    }

    fprintf(file, "%s\n", results);  // Write the results to the file
    fclose(file); // Close the file

    return 0;

}

float calculate_compression_rate_geco (char * sequence_read, int id) {
    FILE *file_seq;
    char filename_uncompressed[100];
    char filename_compressed[100];
    char logs_file[100];
    char command_geco[512];

    // Prepare filenames
    sprintf(filename_uncompressed, "sequence_%d.seq", id);
    sprintf(filename_compressed, "sequence_%d.seq.co", id);
    sprintf(logs_file, "logs_%d.txt", id);

    // Open file for writing
    file_seq = fopen(filename_uncompressed, "w");
    if (!file_seq) {
        perror("Failed to open uncompressed file");
        return -1;
    }

    // Write content
    fprintf(file_seq, "%s\n", sequence_read);
    fclose(file_seq);  // Flush and close before checking size

    // Get size of uncompressed file
    file_seq = fopen(filename_uncompressed, "rb");
    if (!file_seq) {
        perror("Failed to open uncompressed file for size");
        return -1;
    }
    fseek(file_seq, 0, SEEK_END);
    long initial_size = ftell(file_seq);
    fclose(file_seq);

    // Build GeCo3 command
    snprintf(command_geco, sizeof(command_geco),
        "GeCo3 -tm 1:1:0:1:0.9/0:0:0 -tm 7:10:0:1:0/0:0:0 -tm 16:100:1:10:0/3:10:0.9 -lr 0.03 -hs 64 %s > %s 2>&1",
        filename_uncompressed, logs_file);

    // Run the command
    int ret = system(command_geco);
    if (ret != 0) {
        fprintf(stderr, "GeCo3 command failed with return code %d\n", ret);
        return -1;
    }

    // Open compressed file and get size
    file_seq = fopen(filename_compressed, "rb");
    if (!file_seq) {
        perror("Failed to open compressed file");
        return -1;
    }
    fseek(file_seq, 0, SEEK_END);
    long compressed_size = ftell(file_seq);
    fclose(file_seq);

    // Clean up files
    remove(filename_compressed);
    remove(filename_uncompressed);
    remove(logs_file);

    // Return compression ratio

    //printf("%ld  %ld\n\n", compressed_size, initial_size);
    return (float) compressed_size / initial_size;

}

// Tasks to be done by each thread
int worker_task(int index_data_sequence){

    char *content_header;
    char *content_sequence;

    //Initial and end position of the header
    int start_pos_header = data_all_sequences[index_data_sequence].init_header;
    int end_pos_header = data_all_sequences[index_data_sequence].end_header;

    // Get sequence header
    pthread_mutex_lock(&input_file_mutex);
    read_file_partially(start_pos_header, end_pos_header, path_input_file, &content_header);
    char* aux_header = content_header;
    pthread_mutex_unlock(&input_file_mutex);

    // Remove the \n and \t from the headers
    char *read_header = remove_newline_and_tab_characters(aux_header);

    //Initial and end position of the sequence   
    int start_pos_sequence = end_pos_header + 1;  // Start position in the file
    int end_pos_sequence = data_all_sequences[index_data_sequence].end_sequence;   // End position in the file

    // Get sequence
    pthread_mutex_lock(&input_file_mutex);
    read_file_partially(start_pos_sequence, end_pos_sequence, path_input_file, &content_sequence);
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
        results = concatenate_strings(results, float_to_string(size_sequence), 1);
    }

    if (calculate_gc_content == 1){
        float cg_content = (float) data_all_sequences[index_data_sequence].cg_content / data_all_sequences[index_data_sequence].number_bases;
        results = concatenate_strings(results, float_to_string(cg_content), 1);
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
        printf("TODO");
    }
    
    if (compression_geco == 1) {
        float compression_rate_geco = calculate_compression_rate_geco(read_sequence, index_data_sequence);
        results = concatenate_strings(results, float_to_string(compression_rate_geco), 1);

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

    // Read the input file and get relevant info on the sequences
    initial_reading();

    // This may not be necessary
    /*if (number_of_threads > number_sequences){
        number_of_threads = number_sequences;
        printf("Number of threads was set to %d due to there only being %d sequences in the input file.\n", number_of_threads, number_sequences);
    }*/

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

    printf("All threads have finished.\n");
    


}