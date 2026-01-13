#ifndef AUXGENOMECLASS_H
#define AUXGENOMECLASS_H

#include <sys/ioctl.h>
#include <unistd.h>
#include <pthread.h>

typedef struct {
    unsigned long long int id;
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


typedef struct {
    int max_size_seq;
    int number_seqs;
} Info_file;

typedef struct {
    int thread_id;
    FILE *file_output;
} thread_arg_t;

char** append(char* arr[], int n, char* ele);

int get_screen_width();

void progress_bar(int total_tasks);

char *concatenate_strings(char *original, const char *app, int add_tab);

char * int_to_string (int value);

char * float_to_string (float value);

char *remove_newline_and_tab_characters(char *text_to_clean);

int check_if_fa_or_fq (char *file_name, int threads);

Info_file calculate_number_sequences (char* path_input_file);

#endif