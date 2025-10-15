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

#include "genomeclass.h"

/* Minimal globals and declarations so this file can be built as a
    standalone executable (mirrors similar globals in genomeclass.c). */
int help_menu = 0;
int option_index = 0;
FILE *file = NULL;

char *path_input_file = NULL; 
char *output_path = "modified_dataset.fasta";
int length_sequences = 150;
float noise_perc = 0.1;
float mutation_perc = 0.0;
int number_of_threads = 1;
Seq_data *data_all_sequences = NULL;

int number_tasks_assigned = 0;

long long int number_pos_data_sequence = 100000000;
int number_sequences = 0;

pthread_mutex_t input_file_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t output_file_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t task_mutex = PTHREAD_MUTEX_INITIALIZER;

static struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"input", required_argument, 0, 'i'},
    {"output", required_argument, 0, 'o'},
    {"length", required_argument, 0, 'l'},
    {"noise", required_argument, 0, 'n'},
    {"substitutions", required_argument, 0, 's'},
    {"threads", required_argument, 0, 't'},
    {0, 0, 0, 0}
};

// Print help menu
void program_usage(char *prog_path) {
    printf("\nUSAGE: .%s -i <input_fasta> -o <output_fasta> -l <length> -n <noise_percentage> -t <threads>\n\n", strrchr(prog_path, '/'));
    printf("Program options -------------------------------------------------------------------------------------------\n");
    printf("-h, --help\t\tPrints this message\n");
    printf("-i, --input\t\tSet input file (FASTA format).\n");
    printf("-o, --output\t\tSet the output dataset file (FASTA format).\n");
    printf("-l, --length\t\tSet the length of the sequences (default: 150).\n");
    printf("-n, --noise\t\tSet the noise added to the sequences. (default: 0.1).\n");
    printf("-s, --substitutions\tSet the percentage of substitutions to be applied to the sequences (default: 0.0).\n");
    printf("-t, --threads\t\tSets the number of threads.\n");
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
    while ((opt = getopt_long(argc, argv, "hi:o:l:n:s:t:", long_options, &option_index))  != -1) {
        
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
            case 't':
                number_of_threads = atoi(optarg);
                if (number_of_threads < 1) {
                    printf("The argument of option -t must be a positive integer.\n");
                    program_usage(prog_path);
                    return 1;
                }
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
    printf("Lenght of the sequences: %d\n", length_sequences);
    printf("Noise rate: %f\n", noise_perc);
    printf("Mutation rate: %f\n", mutation_perc);
    printf("Number of threads: %d\n", number_of_threads);
 
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

        } else { // Is part of the sequence, increment sequence length and number of bases
            if ((char)ch != '\n') {
                number_bases_seq ++;
            }

            length_seq ++;
        }

        index_file ++;
    }

    // Update the information of the last sequence
    data_all_sequences[index_data_sequences].end_sequence = index_file-1; 
    data_all_sequences[index_data_sequences].length_sequence = length_seq;
    data_all_sequences[index_data_sequences].number_bases = number_bases_seq;

    number_sequences = index_data_sequences + 1;
    printf("Number of sequences in the file - %d\n", number_sequences);   
    
    return 0;

}

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


char *make_changes_sequence(char *sequence_to_mutate){

    char* aux = malloc(strlen(sequence_to_mutate) +1);
    int number_positions = 0;
    int starting_point = (strlen(sequence_to_mutate) - length_sequences);
    int new_start = starting_point * (double)rand() / RAND_MAX;
    int aux_len_sequences = length_sequences;

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

// Tasks to be done by each thread
int worker_task(int index_data_sequence){
   
    char *content_header;
    char *content_sequence;
    char seq_filename[100];

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
    char *out_sequence = make_changes_sequence(aux_sequence);

    // Write to output file
    pthread_mutex_lock(&output_file_mutex);

    FILE *output_file = fopen(output_path, "a");
    if (output_file == NULL) {
        perror("Error opening output file");
        pthread_mutex_unlock(&output_file_mutex);
        return 1;
    }
    fprintf(output_file, "%s\n", read_header);
    fprintf(output_file, "%s\n", out_sequence);
    fclose(output_file);
    pthread_mutex_unlock(&output_file_mutex);
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

    remove(output_path);

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