
#include <sys/ioctl.h>
#include <unistd.h>
#include <pthread.h>
#include <stdio.h>
#include <auxgenomeclass.h>

#define READ_BUF_SIZE (1 << 20)

int tasks_done = 0;
pthread_mutex_t tasks_done_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t progress_bar_mutex = PTHREAD_MUTEX_INITIALIZER;


char** append(char* arr[], int n, char* ele) {
   // Allocate memory for the new array of strings
   char** arrnew = (char**)malloc((n + 1) * sizeof(char*)); 

   // Copy the old strings to the new array
   for (int i = 0; i < n; i++) {
       arrnew[i] = arr[i];  // Copy pointer values (if memory management is correct)
   }

   // Allocate memory for the new string and copy it
   arrnew[n] = (char*)malloc(strlen(ele) + 1);  // +1 for the null terminator
   strcpy(arrnew[n], ele);  // Copy the new string

   return arrnew;
}


int get_screen_width() {
    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    return w.ws_col;  // Use character width, not pixels
}


void progress_bar(int total_tasks) {
    
    pthread_mutex_lock(&tasks_done_mutex);
    tasks_done++;
    pthread_mutex_unlock(&tasks_done_mutex);

    float percentage = (float)tasks_done / total_tasks;
    int percentage_int = (int)(percentage * 100);
    int width = get_screen_width();

    int reserved = 28;
    int bar_width = width - reserved;
    if (bar_width < 10) bar_width = 10; // fallback if terminal is too small

    int pos = (int)(bar_width * percentage);

    pthread_mutex_lock(&progress_bar_mutex);

    printf("\r[");
    for (int i = 0; i < bar_width; i++) {
        if (i < pos) printf("=");
        else if (i == pos) printf(">");
        else printf(" ");
    }
    printf("] %3d%% (%d/%d)", percentage_int, tasks_done, total_tasks);
    fflush(stdout);

    pthread_mutex_unlock(&progress_bar_mutex);
}


char *concatenate_strings(char *original, const char *app, int add_tab) {
    if (!original) original = strdup("");
    if (!app) app = "";

    size_t needed = strlen(original) + strlen(app) + 10;
    char *temp = malloc(needed);

    if (add_tab)
        sprintf(temp, "%s\t%s", original, app);
    else
        sprintf(temp, "%s%s", original, app);

    free(original);  // THIS IS CRITICAL

    return temp;
}


char * int_to_string (int value) {

    char *temp = malloc(10);
    sprintf(temp, "%d", value);   
    return temp;

}


char * float_to_string (float value) {

    char *temp = malloc(32);
    sprintf(temp, "%f", value);   
    return temp;

}


char *remove_newline_and_tab_characters(char *text_to_clean){

    char* aux = malloc(strlen(text_to_clean) +1);
    int number_positions = 0;

    for (size_t i = 0; i < strlen(text_to_clean); i++) {
        if (text_to_clean[i] != '\n' && text_to_clean[i] != '\t'){ // Remove \n and \t from the headers
            aux[number_positions] = text_to_clean[i];
            number_positions++;
        }
    }
    aux[number_positions] = '\0';

    return aux;

}


int check_if_fa_or_fq (char *file_name, int threads) {

    char command_spades[512];
    char new_path[100];
    char *dot = strrchr(file_name, '.');  // find last occurrence of '.'

    if (dot != NULL && *(dot + 1) != '\0') {
        //printf("\n\nAfter last . :-%s-\n", dot + 1);

        if (strcmp(dot+1, "fasta") == 0 || strcmp(dot+1, "fa") == 0) {
            printf("FASTA file detected\n");
            return 0;

        } else if (strcmp(dot+1, "fastq") == 0 || strcmp(dot+1, "fq") == 0) {
            printf("FASTQ file detected - reconstructing the sample\n");

            
            // Build spades command 
            snprintf(command_spades, sizeof(command_spades), "conda run -n genomeclass spades.py -s %s -o spades_output -t %d", file_name, threads);

            
            // Run the spades command
            int ret = system(command_spades);
            if (ret != 0) {
                fprintf(stderr, "spades command failed with return code %d\n", ret);
                return -1;
            }

            if (fopen("spades_output/contigs.fasta", "r")) {
                printf("Contig file exists!\n");
                snprintf(new_path, sizeof(new_path), "mv spades_output/contigs.fasta %s", file_name);
                system(new_path);
                remove("spades_output");
                return 0;
            } else if (fopen("spades_output/scaffolds.fasta", "r")) {
                printf("File does not exist.\n");
                snprintf(new_path, sizeof(new_path), "mv spades_output/scaffolds.fasta %s", file_name);
                system(new_path);
                remove("spades_output");
                return 0;
                
            } else {
                return 1;
            }

        } else {
            printf("Error - The input file is not supported.\n");
            return 1;
        }
    } else {
        printf("Error - The input file is not supported.\n");
        return 1;
    }

}


Info_file calculate_number_sequences(char *path_input_file) {

    Info_file info_file;
    memset(&info_file, 0, sizeof(info_file));

    FILE *file = fopen(path_input_file, "r");
    if (!file) {
        perror("Error opening file");
        return info_file;
    }

    char buffer[READ_BUF_SIZE];
    size_t bytes_read;

    int is_header = 0;
    size_t count_bases = 0;

    while ((bytes_read = fread(buffer, 1, READ_BUF_SIZE, file)) > 0) {

        for (size_t i = 0; i < bytes_read; ++i) {
            char ch = buffer[i];

            if (is_header) {
                if (ch == '\n') {
                    is_header = 0;
                }
            }
            else {
                if (ch == '>') {
                    info_file.number_seqs++;

                    if (count_bases > info_file.max_size_seq) {
                        info_file.max_size_seq = count_bases;
                    }

                    count_bases = 0;
                    is_header = 1;
                }
                else if (ch != '\n') {
                    count_bases++;
                }
            }
        }
    }

    /* Final sequence (EOF without trailing '>') */
    if (count_bases > info_file.max_size_seq) {
        info_file.max_size_seq = count_bases;
    }

    fclose(file);

    return info_file;
}
