#include <htslib/sam.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstring>

using namespace std;

char** chr_table = NULL;
int n_chr = 0;

// only use malloc()/free() for chromosome table
char* get_chr_ptr(const char* chr_name){
	if(!chr_table) chr_table = (char**)malloc(sizeof(char*)); 
	for(int i = 0; i < n_chr; ++i){
		if(!strcmp(chr_name, chr_table[i])){
			return chr_table[i];
		}
	}
	// didn't find it in chr_table, so allocate another 
	chr_table = (char**)realloc(chr_table, sizeof(char*)*(++n_chr));
	chr_table[n_chr-1] = (char*)malloc(sizeof(char)*strlen(chr_name));
	strcpy(chr_table[n_chr-1], chr_name);
	return chr_table[n_chr-1];
}

void free_chr_table(){
	for(int i = 0; i < n_chr; ++i) free(chr_table[i]);
	free(chr_table);
}

typedef struct read_t{
	int id;
	char* chr;
	int pos;
	int strand;
	int ins;
}read_t;

typedef struct read_range{
	read_t* begin;
	read_t* end;
}read_range;

read_t* read_arr;
read_range* peak_idx;

int main(int argc, char* argv[]){
	/**
	 * Usage: main [reads] [no. peaks] [no. reads]
	 */

	const char* reads_file_name = argv[1];
	int n_peaks = atoi(argv[2]);
	int n_reads = atoi(argv[3]);

	read_arr = new read_t[n_reads];
	peak_idx = new read_range[n_peaks];

	ifstream reads_file(reads_file_name);
	string dummy; std::getline(reads_file, dummy);
	
	int index = 0;
	int current_peak_id = -1;
	
	while(!reads_file.eof()){
		
		string id, chr, pos, strand, ins;
		if(std::getline(reads_file, id, '\t') &&\
			std::getline(reads_file, chr, '\t') &&\
			std::getline(reads_file, pos, '\t') &&\
			std::getline(reads_file, strand, '\t') &&\
			std::getline(reads_file, ins, '\n')){
		}else{
			break;
		}
		int id_int = stoi(id);
		int pos_int = stoi(pos);
		int strand_int = strand[0] == '+' ? 1 : -1;
		int ins_int = stoi(ins);
		read_arr[index] = (read_t){.id = id_int, .chr = get_chr_ptr(chr.c_str()), .pos=pos_int, .strand=strand_int, .ins = ins_int};
		if(current_peak_id == -1){
			current_peak_id = id_int;
			peak_idx[0].begin = &read_arr[index];
		}else if(current_peak_id != id_int){
			peak_idx[current_peak_id-1].end = index > 0 ? &read_arr[index-1] : &read_arr[0];
			peak_idx[current_peak_id].begin = &read_arr[index];
			current_peak_id = id_int;
		}
		index++;
		if(current_peak_id > 10) break;
	}
	free_chr_table();
	return 0;
}
