#include <htslib/sam.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstring>
#include <pthread.h>
#include <vector>
#include <unistd.h>
#include <climits>
#include <algorithm>

#include "kde/kde.h"
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

/** using some global variables **/

struct {
	double bandwidth = 12;
	double abs_thresh = 0;
	double rel_thresh = 0.5;
	double prominence = 0.9;
}param;

double density(read_range reads, double bw, double* x, double* y, int n=512){
	// returns the max y-value
	// because we need it later!
	read_t* r = reads.begin;
	int n_reads = (reads.end - reads.begin) + 1;
	/** bunch of internal parameters (might expose later) **/
	double cut = 3;
	int min = INT_MAX; int max = INT_MIN;

	do{
		if(r->pos > max) max = r->pos;
		if(r->pos < min) min = r->pos;
	}while(r++ != reads.end); /** determine the range of pos values **/

	/** behave like density() from R **/

	double start = min - cut*bw;
	double end = max + cut*bw; // see help for density() in R
	double step = (end - start)/n;

	double max_y = 0.0;

	int i = 0;
	for(double obs = start; obs < end - (step/2); obs += step){
		double prob = 0;
		r = reads.begin;
		do{
			/** calculate kernel density **/
			prob += gauss_kernel( (r->pos - obs)/bw)/(n_reads*bw);
		}while(r++ != reads.end);
		x[i] = obs;
		y[i] = prob;
		if(prob > max_y){
			max_y = prob;
		}
		i++;
	}
	return max_y;
}

vector<int> find_peaks(double* y, int start, int end, double thresh){
	/** find_peaks will work on [start, end)
	**/
	
	if(start >= end) return vector<int>();

	double max = -1; int max_idx = 0;
	for(int i = start; i<end; i++){
		if(y[i] > max){
			max = y[i]; max_idx = i;
		}
	}
	
	if(max < thresh) return vector<int>();

	vector<int> peaks = {max_idx};

	/** descend forward from the scope global max **/
	int i = max_idx;
	while(i < end && y[i] > max*param.prominence) i++;
	double min = y[i];
	while(i < end && y[i] < min/param.prominence){
		if(min > y[i]) min=y[i]; i++;
	}
	/** now RECURSE :o **/
	vector<int> peaks2 = find_peaks(y, i, end, thresh);
	peaks.insert(peaks.end(), peaks2.begin(), peaks2.end());

	/** descend backward from the scope global max **/
	i = max_idx;
	while(i >= start && y[i] > max*param.prominence) i--;
	min = y[i];
	while(i >= start && y[i] < min/param.prominence){
		if(min > y[i]) min = y[i]; i--;
	}
	/** now RECURSE again OwO **/
	vector<int> peaks3 = find_peaks(y, start, i, thresh);
	peaks.insert(peaks.begin(), peaks3.begin(), peaks3.end());

	return peaks;
}

void* thread_worker(void* arg){
	/** thread worker **/
	/** arg --> std::pair<int, int> that 
	 * specify the first and last peaks that contribute
	 * to given workload
	 */
	
	int a, b;
	a = ((std::pair<int, int>*)arg)->first;
	b = ((std::pair<int, int>*)arg)->second;


	double x[512];
	double y[512];

	for(int i = a; i <= b; i++){
		// insert - check for n_reads > threshold
		double max_y = density(peak_idx[i], param.bandwidth, x, y);
		double threshold = std::max(param.abs_thresh, max_y*param.rel_thresh);
		vector<int> peaks = find_peaks(y, 0, 512, threshold);
		cout << "size = " << peaks.size() << endl;
	}
		
	return NULL;
}

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

	/** code for multithreading **/

	int n_threads = 1;

	pthread_t* threads;
	void** thread_args;
	void** thread_retvals;

	threads = new pthread_t[n_threads];
	thread_retvals = new void*[n_threads];
	thread_args = new void*[n_threads];

	
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
			peak_idx[current_peak_id-1].begin = &read_arr[index];
		}else if(current_peak_id != id_int){
			current_peak_id = id_int;
			peak_idx[current_peak_id-2].end = index > 0 ? &read_arr[index-1] : &read_arr[0];
			peak_idx[current_peak_id-1].begin = &read_arr[index];
		}
		index++;
	}
	peak_idx[current_peak_id-1].end = index > 0 ? &read_arr[index-1] : &read_arr[0];

	/** assign peaks **/
	/* 
	 *  thread_args will contain a std::pair<int, int> corresponding to [start_peak, end_peak]
	 */

	for(int i = 0; i < n_peaks; i++){
		cout << "Size: " << peak_idx[i].end - peak_idx[i].begin << endl;
	}
	
	int split_size = n_peaks/n_threads;
	int start = 0;	
	for(int i = 0; i < n_threads; i++){
		int end = (i == n_threads - 1) ? n_peaks-1 : start + split_size-1;
		thread_args[i] = new std::pair<int, int>(start, end); // this is an inclusive range
		start = end+1;
	}


	/** spawn worker threads **/
	for(int i = 0; i < n_threads; i++){
		cout << "Creating thread " << i << endl;
		pthread_create(&threads[i], NULL, thread_worker, thread_args[i]);
	}

	for(int i = 0; i < n_threads; i++){
		pthread_join(threads[i], NULL);
	}

	delete[] threads;
	
	for(int i = 0; i < n_threads; i++){
		delete (std::pair<int, int>*)thread_args[i];
	}

	delete[] thread_args;
	delete[] thread_retvals;
	
	free_chr_table();
	return 0;
}
