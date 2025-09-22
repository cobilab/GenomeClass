
#include <sys/ioctl.h>
#include <unistd.h>
// Structure for storing the beginning and end of a sequence

typedef struct {
    unsigned long long int init_header;
    unsigned long long int end_header;
    unsigned long long int end_sequence;
    unsigned long long int length_sequence;
    unsigned long long int number_bases;
    unsigned long long int cg_content;
    unsigned long long int number_a;
    unsigned long long int number_c;
    unsigned long long int number_g;
    unsigned long long int number_t;
    unsigned long long int number_other;
} Seq_data;

typedef struct {
    float avg_distance;
    float prob_sequence;
} Dist_Prob_sequence;

int tasks_done = 0;
pthread_mutex_t tasks_done_mutex = PTHREAD_MUTEX_INITIALIZER;


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

    float percentage = (float)tasks_done / total_tasks;
    int percentage_int = (int)(percentage * 100);
    int width = get_screen_width();

    // Reserve space for " 100% (xxxx/xxxx)" ≈ 20 chars
    int reserved = 20;
    int bar_width = width - reserved;
    if (bar_width < 10) bar_width = 10; // fallback if terminal is too small

    int pos = (int)(bar_width * percentage);

    printf("\r[");
    for (int i = 0; i < bar_width; i++) {
        if (i < pos) printf("=");
        else if (i == pos) printf(">");
        else printf(" ");
    }
    printf("] %3d%% (%d/%d)", percentage_int, tasks_done, total_tasks);
    fflush(stdout);

    pthread_mutex_unlock(&tasks_done_mutex);
}



char * concatenate_strings(char *original_content, char *content_to_append, int add_tab) {

    char *temp = malloc(strlen(original_content) + strlen(content_to_append) + 10);
    if (add_tab == 0){
        sprintf(temp, "%s%s", original_content, content_to_append);
    } else {
        sprintf(temp, "%s\t%s", original_content, content_to_append);
    }
    
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


