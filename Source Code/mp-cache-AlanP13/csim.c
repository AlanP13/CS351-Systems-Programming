/******************************************************************************
 *
 * Cache lab header
 *
 * Author: Alan Palayil
 * Email:  apalayil@hawk.iit.edu
 * AID:    A20447935
 * Date:   10/25/21
 *
 * By signing above, I pledge on my honor that I neither gave nor received any
 * unauthorized assistance on the code contained in this repository.
 *
 *****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>
#include "cachelab.h" 
const int m = sizeof(long) * 8; 
int hit_count = 0, miss_count = 0, evict_count = 0;
int v_flag = 0, h_flag = 0;
int s, b, E, t;
char *trace;
struct Line {
  int valid ;
  unsigned tag ;
  struct Line *next;
};

struct Set {
  int E;
  struct Line *head_line;
  struct Set *next;
};

struct Cache {
  int S;
  struct Set *head_set;
};
int set_size(struct Set *set){
  struct Line *curr = set->head_line;
  int size = 0;
  while (curr != NULL){
    ++size;
    curr = curr->next;
  }
  return size;
}
void init_set(struct Set *set, int E){
  set->E = E;
  set->head_line = NULL;
  set->next = NULL;
}
void init_cache(struct Cache *cache, int s, int E) {
  cache->S = (1 << s);
  struct Set *thisset = (struct Set*)malloc(sizeof(struct Set));
  init_set(thisset, E);
  cache->head_set= thisset; 
  int i;
  for (i = 1; i < cache->S; ++i){
    struct Set *nextset = (struct Set*)malloc(sizeof(struct Set));
    init_set(nextset, E);
    thisset->next = nextset; 
    thisset = nextset;
  }
}
void delete_last_line(struct Set *set){
  struct Line *curr = set->head_line;
  struct Line *prev = NULL;
  struct Line *prev2 = NULL;
  while (curr != NULL){
    prev2 = prev;
    prev = curr;
    curr = curr->next;
  }
  if (prev2 != NULL){
    prev2->next = NULL;
  } else {
    set->head_line = NULL;
  }
  if (prev != NULL) free(prev);
}
void add_line_to_head(struct Set *set, struct Line *line){
  if (set_size(set) == set->E){
    delete_last_line(set);
    evict_count++;
    if (v_flag) printf("eviction ");
  }
  line->next = set->head_line;
  set->head_line = line;
}
void move_line_to_head(struct Set *set, struct Line *line, struct Line *prev){
  if (prev != NULL){
    prev->next = line->next;
    line->next = set->head_line;
    set->head_line = line;
  }
}
void access_cache(struct Cache *cache, unsigned address){
  int tag_bits = address >> (s+b);
  int set_bits = (address << t) >> (t+b);
  struct Set *target_set = cache->head_set;
  int i;
  for (i = 0; i < set_bits; ++i){
    target_set = target_set->next;
  } 
  struct Line *ln = target_set->head_line;
  struct Line *prev = NULL;
  while (ln != NULL){
     if (ln->valid && (ln->tag == tag_bits)){
      hit_count++;
      if (v_flag) printf("hit ");
      move_line_to_head(target_set, ln, prev);
      return ;
    }
    prev = ln;
    ln = ln->next;
  }
  miss_count++;
  if (v_flag) printf("miss ");
  struct Line *newln = (struct Line *)malloc(sizeof(struct Line));
  newln->valid = 1;
  newln->tag = tag_bits;
  add_line_to_head(target_set,newln);
}
void free_cache(struct Cache *cache){
  struct Set *set_to_free = cache->head_set;
  while (!set_to_free){
    struct Line *ln_to_free = set_to_free->head_line;
    while (!ln_to_free){
      struct Line *temp_ln = ln_to_free->next;
      free(ln_to_free);
      ln_to_free = temp_ln;
    }
    struct Set *temp_set = set_to_free->next;
    free(set_to_free);
    set_to_free = temp_set;
  }
  free(cache);
}
int main(int argc, char** argv){
  int opt;
  char op;
  unsigned addr;
  int size;
  while(-1 != (opt = getopt(argc, argv, "vhs:E:b:t:"))) {
    switch(opt) {
    case 'v':
      v_flag = 1;
      break;
    case 'h':
      h_flag = 1;
      break;
    case 's':
      s = atoi(optarg);
      if(s < 0 && s > m ){
	printf("Error: s value must be in [0, word_length].\n");
	exit(-1);
      }
      break;
    case 'E':
      E = atoi(optarg);
      if(E <= 0){
	printf("Error: E value must be larger than 0.\n");
	exit(-1);
      }
      break;
    case 'b':
      b = atoi(optarg);
      if(b < 0 || b > m){
	printf("Error: b value must be in [0, word_length].\n");
	exit(-1);
      }
      break;
    case 't':
      trace = optarg;
      break;
    default:
      break;
    }
  }
  t = m - s - b;
  struct Cache *mycache = (struct Cache *)malloc(sizeof(struct Cache));
  init_cache(mycache, s, E);
  FILE *traceFile;
  traceFile = fopen(trace, "r");
  if (!traceFile){
    fprintf(stderr, "Error: Trace file cannot be opened.\n");
    return -1;
  } 
  while(fscanf(traceFile, "%c %x, %d", &op, &addr, &size) > 0) {
    if (v_flag) printf("%c %x, %d ", op, addr, size);
    switch(op){
    case 'L':
      access_cache(mycache, addr);
      break;
    case 'S':
      access_cache(mycache, addr);
      break;
    case 'M':
      access_cache(mycache, addr);
      access_cache(mycache, addr);
      break;
    }
    if (v_flag) printf("\n");
  }
  fclose(traceFile); 

  free_cache(mycache);
  printSummary(hit_count, miss_count, evict_count);

  return 0; 
}
