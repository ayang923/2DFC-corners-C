#ifndef __HASHMAP_H__
#define __HASHMAP_H__

#include <stdio.h>
#include <stdbool.h>
#include <string.h>

typedef struct bucket{
    int key;
    double value;
    bool is_occupied;
} bucket_t;

typedef struct hashmap{
    bucket_t *buckets; // Fixed-size bucket array
    size_t max_buckets;
} hashmap_t;

void hashmap_init(hashmap_t *map, size_t max_buckets);

bool hashmap_insert(hashmap_t *map, int key, double value);

bool hashmap_get(hashmap_t *map, int key, double *value);

#endif