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
#include <sys/syscall.h>
#include <sys/types.h>
#include <argparse/argparse.hpp>

#include "kde/kde.h"

using namespace std;

char** chr_table = NULL;
int n_chr = 0;

typedef struct footprint_t{
	vector<int> peaks;
	vector<int> fp;
	int centre;
	int peakidx;
}footprint_t;

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
	int min_reads = 500;
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
	for(double obs = start; i < n; obs += step){
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
	while(i > start && y[i] > max*param.prominence) i--;
	min = y[i];
	while(i > start && y[i] < min/param.prominence){
		if(min > y[i]) min = y[i]; i--;
	}
	/** now RECURSE again OwO **/
	vector<int> peaks3 = find_peaks(y, start, i, thresh);
	peaks.insert(peaks.begin(), peaks3.begin(), peaks3.end());

	return peaks;
}

vector<int> find_footprints(double* y, vector<int>& peaks){
	if(peaks.size() > 1){
		vector<int> fp(peaks.size()-1);
		for(int i = 0; i < peaks.size()-1; i++){
			double min = DBL_MAX; int min_idx = peaks[i];
			for(int j = peaks[i]; j < peaks[i+1]; j++){
				if(y[j] < min){
					min = y[j]; min_idx = j; 
				}
			}
			fp[i] = min_idx; 
		}
		return fp;
	}
	return vector<int>();	// no footprints
}

int find_centre(vector<int>& peaks, vector<int>& fp){
	if(fp.size() % 2 == 0){ // even # of footprints, centre on middle peak
		return peaks[peaks.size()/2];
	}
	// odd # of footprints, centre on middle footprint
	return fp[fp.size()/2];
}

void* thread_worker(void* arg){
	/** thread worker **/
	/** arg --> std::pair<int, int> that 
	 * specify the first and last peaks that contribute
	 * to given workload
	 */
	cerr << "Thread starting: " << (long int)syscall(SYS_gettid) << endl;
	
	int a, b;
	a = ((std::pair<int, int>*)arg)->first;
	b = ((std::pair<int, int>*)arg)->second;

	vector<footprint_t>* results = new vector<footprint_t>; results->reserve(b - a + 1);

	double x[512];
	double y[512];

	for(int i = a; i <= b; i++){
		int n_reads = peak_idx[i].end - peak_idx[i].begin + 1;
		if(n_reads < param.min_reads) continue;
		
		double max_y = density(peak_idx[i], param.bandwidth, x, y);
		double threshold = std::max(param.abs_thresh, max_y*param.rel_thresh);
		vector<int> peaks = find_peaks(y, 0, 512, threshold);
		if(peaks.size() > 0){
			/* get footprint positions */
			vector<int> fp = find_footprints(y, peaks);
			if(fp.size() > 0){
				vector<int> fp_pos(fp.size());
				for(int i = 0; i < fp.size(); i++){
					fp_pos[i] = (int)x[fp[i]];
				}
				// centres
				int centre = (int)x[find_centre(peaks, fp)];
				results->push_back((footprint_t){.peaks = peaks, .fp = fp_pos, .centre = centre, .peakidx = i}); 
			}
		}else{

		}
	}
	cerr << "Thread terminating: "<< syscall(SYS_gettid) << endl;		
	return (void*)results;
}

int main(int argc, const char* argv[]){
	/**
	 * Usage: main [reads] [no. peaks] [no. reads]
	 */
	ArgumentParser parser;
	parser.addArgument("-N", "--nreads", 1);
	parser.addArgument("-P", "--npeaks", 1);
	parser.addArgument("-B", "--bw", 1);
	parser.addArgument("-T", "--absthresh", 1);
	parser.addArgument("-t", "--relthresh", 1);
	parser.addArgument("-p", "--prominence", 1);
	parser.addArgument("-m", "--minreads", 1);
	parser.addArgument("-n", "--threads", 1);
	parser.addFinalArgument("input");
	
	parser.parse(argc, argv);
	std::string reads_file_name_str = parser.retrieve<string>("input");
	const char* reads_file_name = reads_file_name_str.c_str();
	
	int n_reads = stoi(parser.retrieve<string>("nreads"));
	int n_peaks = stoi(parser.retrieve<string>("npeaks"));

	if(parser.retrieve<string>("bw") != "") param.bandwidth = stof(parser.retrieve<string>("bw"));
	if(parser.retrieve<string>("absthresh") != "" ) param.abs_thresh = stof(parser.retrieve<string>("absthresh"));
	if(parser.retrieve<string>("relthresh") != "" ) param.rel_thresh = stof(parser.retrieve<string>("relthresh"));
	if(parser.retrieve<string>("prominence") != "" ) param.prominence = stof(parser.retrieve<string>("prominence"));
	if(parser.retrieve<string>("minreads") != "" ) param.min_reads = stof(parser.retrieve<string>("minreads"));
	
	read_arr = new read_t[n_reads];
	peak_idx = new read_range[n_peaks];

	ifstream reads_file(reads_file_name);
	string dummy; std::getline(reads_file, dummy);
	
	int index = 0;
	int current_peak_id = -1;

	/** code for multithreading **/

	int n_threads = 1;
	if(parser.retrieve<string>("threads") != "") n_threads = stoi(parser.retrieve<string>("threads")); 

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

	// for(int i = 0; i < n_peaks; i++){
	// 	cout << "Size: " << peak_idx[i].end - peak_idx[i].begin << endl;
	// }
	
	int split_size = n_peaks/n_threads;
	int start = 0;	
	for(int i = 0; i < n_threads; i++){
		int end = (i == n_threads - 1) ? n_peaks-1 : start + split_size-1;
		thread_args[i] = new std::pair<int, int>(start, end); // this is an inclusive range
		start = end+1;
	}


	/** spawn worker threads **/
	for(int i = 0; i < n_threads; i++){
		pthread_create(&threads[i], NULL, thread_worker, thread_args[i]);
	}

	for(int i = 0; i < n_threads; i++){
		pthread_join(threads[i], &thread_retvals[i]);
	}

	cerr << "Finished processing peaks" << endl;

	cout << "peakid n_peaks n_fp pos centre" << endl;
	for(int i = 0; i < n_threads; i++){
		// do somethign with output...
		vector<footprint_t> *res = (vector<footprint_t>*)thread_retvals[i];
		for(auto i:*res){
			for(auto j:i.fp){
				cout << i.peakidx+1 << "\t" << i.peaks.size() << "\t" << i.fp.size() << "\t" << j << "\t" << i.centre << endl;
			}
		}
		delete res;
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
