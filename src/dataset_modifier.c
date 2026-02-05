#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <time.h>
#include <ctype.h>
#include <pthread.h>

#include "auxgenomeclass.h"

/* Minimal globals and declarations so this file can be built as a
    standalone executable (mirrors similar globals in genomeclass.c). */
int help_menu = 0;
int option_index = 0;
FILE *file = NULL;

char *path_input_file = NULL;
char *output_path = "modified_dataset.fasta";
int length_sequences = -1;
float noise_perc = 0.1;
float mutation_perc = 0.0;
int number_of_threads = 1;
Seq_data *data_all_sequences = NULL;

int number_tasks_assigned = 0;

long long int number_pos_data_sequence = 100000000;
int number_sequences = 0;
int max_number_bases = 0;
int verbose = 0;

float max_number_sequences = 1.0;

pthread_mutex_t input_file_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t output_file_mutex = PTHREAD_MUTEX_INITIALIZER;

static struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"input", required_argument, 0, 'i'},
    {"output", required_argument, 0, 'o'},
    {"length", required_argument, 0, 'l'},
    {"noise", required_argument, 0, 'n'},
    {"substitutions", required_argument, 0, 's'},
    {"max_sequences", required_argument, 0, 'm'},
    {"threads", required_argument, 0, 't'},
    {"disable_pbar", no_argument, 0, 'd'},
    {0, 0, 0, 0}
};

// Print help menu
void program_usage(char *prog_path) {
    printf("\nUSAGE: .%s -i <input_fasta> -o <output_fasta> -l <length> -n <noise_percentage> -t <threads>\n\n", strrchr(prog_path, '/'));
    printf("Program options -------------------------------------------------------------------------------------------\n");
    printf("-h, --help\t\tPrints this message\n");
    printf("-i, --input\t\tSet input file (FASTA format).\n");
    printf("-o, --output\t\tSet the output dataset file (FASTA format).\n");
    printf("-l, --length\t\tSet the maximum length of the sequences (default: not set).\n");
    printf("-n, --noise\t\tSet the noise added to the sequences. (default: 0.1).\n");
    printf("-s, --substitutions\tSet the percentage of substitutions to be applied to the sequences (default: 0.0).\n");
    printf("-m, --max_sequences\tSet the percentage of sequences to be output.\n");
    printf("-t, --threads\t\tSets the number of threads.\n");
    printf("-d, --disable_pbar\tDisable progress bar output.\n");
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
    while ((opt = getopt_long(argc, argv, "hi:o:l:n:s:m:t:d", long_options, &option_index))  != -1) {

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
            case 'l':
                length_sequences = atoi(optarg);
                if (length_sequences < 1) {
                    printf("The argument of option -l must be a positive integer.\n");
                    program_usage(prog_path);
                    return 1;
                }
                break;
            case 'n':
                noise_perc = atof(optarg);
                if (noise_perc < 0 || noise_perc > 1) {
                    printf("The argument of option -n must be a float between 0 and 1.\n");
                    program_usage(prog_path);
                    return 1;
                }
                break;
            case 's':
                mutation_perc = atof(optarg);
                if (mutation_perc < 0 || mutation_perc > 1) {
                    printf("The argument of option -m must be a float between 0 and 1.\n");
                    program_usage(prog_path);
                    return 1;
                }
                break;
            case 'm':
                max_number_sequences = atof(optarg);
                if (max_number_sequences < 0 || max_number_sequences > 1) {
                    printf("The argument of option -m must be a float between 0 and 1.\n");
                    program_usage(prog_path);
                    return 1;
                }
                break;
            case 't':
                number_of_threads = atoi(optarg);
                if (number_of_threads < 1) {
                    printf("The argument of option -t must be a positive integer.\n");
                    program_usage(prog_path);
                    return 1;
                }
                break;
            case 'd':
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
    printf("Input file: %s\n", path_input_file ? path_input_file : "None");
    printf("Output file: %s\n", output_path ? output_path : "None");
    printf("Maximum length of the sequences: %d\n", length_sequences);
    printf("Noise rate: %f\n", noise_perc);
    printf("Mutation rate: %f\n", mutation_perc);
    printf("Ratio of sequences selected: %f\n", max_number_sequences);
    printf("Number of threads: %d\n", number_of_threads);

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

char *make_changes_sequence(char *sequence_to_mutate){

    char* aux = malloc(strlen(sequence_to_mutate) +1);
    int number_positions = 0;

    int starting_point = 0;
    int new_start = 0;
    int aux_len_sequences = strlen(sequence_to_mutate);

    if (length_sequences > 0 && strlen(sequence_to_mutate) >= length_sequences) {
        starting_point = (strlen(sequence_to_mutate) - length_sequences);
        new_start = starting_point * (double)rand() / RAND_MAX;
        aux_len_sequences = length_sequences;
    }

    for (size_t i = 0; i < strlen(sequence_to_mutate); i++) {

        if (i < new_start || i >= new_start + aux_len_sequences){ // Skip the bases outside the desired length
            continue;
        } else {

            if (sequence_to_mutate[i] != '\n' && sequence_to_mutate[i] != '\t'){ // Remove \n and \t from the headers

                if ((char)sequence_to_mutate[i] == 'a' || (char)sequence_to_mutate[i] == 'c' || (char)sequence_to_mutate[i] == 't' || (char)sequence_to_mutate[i] == 'g'|| (char)sequence_to_mutate[i] == 'n'){
                    sequence_to_mutate[i] = (char) toupper(sequence_to_mutate[i]);
                }

                if (sequence_to_mutate[i] == 'N' ){
                    aux[number_positions] = sequence_to_mutate[i];
                } else if ((float)rand() / RAND_MAX < noise_perc){ // If the random value is less than the noise percentage, mutate the base
                    aux[number_positions] = 'N';

                } else if ((float)rand() / RAND_MAX < mutation_perc){ // If the random value is less than the mutation percentage, mutate the base
                    const char valid_bases[] = "ACTG";
                    aux[number_positions] = valid_bases[rand() % 4];
                } else {
                    aux[number_positions] = sequence_to_mutate[i];
                }

                number_positions++;
            } else {
                aux_len_sequences++;
            }
        }
    }
    aux[number_positions] = '\0';

    return aux;

}



char* make_sequence(unsigned long long int start_pos_sequence, unsigned long long int end_pos_sequence) {

    char *content_sequence;
    char *processed_sequence;

    // Get sequence
    read_file_partially(start_pos_sequence, end_pos_sequence, &content_sequence);
    char* aux_sequence = content_sequence;

    processed_sequence = make_changes_sequence(aux_sequence);

    free(aux_sequence);

    return processed_sequence;
}

int dump_sequence(char* header, char* results) {

    if (results != NULL) {

        // Write to output file
        pthread_mutex_lock(&output_file_mutex);

        FILE *output_file = fopen(output_path, "a");
        if (output_file == NULL) {
            perror("Error opening output file");
            pthread_mutex_unlock(&output_file_mutex);
            return 1;
        }
        fprintf(output_file, "%s\n", header);
        fprintf(output_file, "%s\n", results);
        fclose(output_file);
        pthread_mutex_unlock(&output_file_mutex);

        free(results);
    }

    return 0;
}


//Process the FASTA file and write output to TSV file
//TODO - check this function and change its logic to work with the new way of processing sequences
int process_file (int thread_id) {

    int ch;
    long int count_sequences = 0;
    int processing = 0;
    int is_header = 0;
    int is_sequence_processed = 0;

    int max_header = 500;
    char header[max_header];
    int pos_header = 0;

    unsigned long long int start_pos_sequence = 0;
    unsigned long long int end_pos_sequence = 0;

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

                if ((float)rand() / RAND_MAX < max_number_sequences) { // Process the content based on the max number of sequences

                    char* results = make_sequence(start_pos_sequence, end_pos_sequence);
                    dump_sequence(header, results);

                }

                // Update progress bar
                if (verbose == 0){ // Update the progress bar
                    progress_bar(number_sequences);
                }
                
                pos_header = 0;
            }
            
            
            count_sequences ++;
            is_header = 1;

            if (count_sequences % number_of_threads == thread_id - 1) { //Process the content
                processing = 1;
                pos_header = 0;
                memset(header, 0, max_header);
                header[0] = '>';
                pos_header ++;
                header[pos_header] = '\0';

            } else {
                processing = 0;
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

                if (is_sequence_processed == 0) {
                    // Update the end position of the sequence
                    end_pos_sequence = current_pos;
                }
                
            }

        }

        current_pos ++;

    }

    // Process the last sequence in the file
    if (processing == 1) { // Was processing info; dump onto file
        
        if ((float)rand() / RAND_MAX < max_number_sequences) { // Process the content based on the max number of sequences

            char* results = make_sequence(start_pos_sequence, end_pos_sequence);
            dump_sequence(header, results);

        }

        // Update progress bar
        if (verbose == 0){ // Update the progress bar
            progress_bar(number_sequences);
        }
    }


    //printf("\nThread %d has finished processing.\n", thread_id);

    fclose(file);

    return 0;
}




// Thread function
void* thread_function(void* arg) {
    
    int thread_number = *(int*)arg + 1;

    process_file(thread_number);

    return NULL;
}





int main(int argc, char *argv[]) {

    srand(time(NULL));  // Seed the random number generator

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

    if (verbose == 0) {
    
        Info_file info_file = calculate_number_sequences(path_input_file);

        number_sequences = info_file.number_seqs;
        max_number_bases = info_file.max_size_seq;
    }

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

}
