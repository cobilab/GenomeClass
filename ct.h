
#include <sys/ioctl.h>
#include <unistd.h>
// Structure for storing the begining and end of a sequence

typedef struct {
    int init_header;
    int end_header;
    int end_sequence;
    int length_sequence;
    int number_bases;
    int cg_content;
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

    int pos = (int)(width * percentage);

    printf("\r[");
    for (int i = 0; i < width - 7; i++) {  // -7 to make room for " 100%"
        if (i < pos) printf("=");
        else if (i == pos) printf(">");
        else printf(" ");
    }
    printf("] %3d%%", percentage_int);
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

    for (int i = 0; i < strlen(text_to_clean); i++){
        if (text_to_clean[i] != '\n' && text_to_clean[i] != '\t'){ // Remove \n and \t from the headers
            aux[number_positions] = text_to_clean[i];
            number_positions++;
        }
    }
    aux[number_positions] = '\0';

    return aux;

}


