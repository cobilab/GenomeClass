
#include <sys/ioctl.h>
#include <unistd.h>
// Structure for storing the begining and end of a sequence

typedef struct {
    int init_header;
    int end_header;
    int end_sequence;
    int length_sequence;
    int number_bases;
} Seq_data;

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

