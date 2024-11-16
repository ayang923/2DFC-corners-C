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

// Simple hash function
size_t hash_function(hashmap_t *map, int key) {
    return key % map->max_buckets;
}

void hashmap_init(hashmap_t *map, size_t max_buckets) {
    memset((int*) map->buckets, 0, max_buckets*sizeof(bucket_t));
    map->max_buckets = max_buckets;
}

// Insert key-value pair into the hashmap
bool hashmap_insert(hashmap_t *map, int key, double value) {
    size_t index = hash_function(map, key);
    for (size_t i = 0; i < map->max_buckets; i++) {
        size_t probe_index = (index + i) % map->max_buckets;
        if (!map->buckets[probe_index].is_occupied || map->buckets[probe_index].key == key) {
            map->buckets[probe_index].key = key;
            map->buckets[probe_index].value = value;
            map->buckets[probe_index].is_occupied = true;
            return true;
        }
    }
    return false; // Hashmap is full
}

// Retrieve a value by key
bool hashmap_get(hashmap_t *map, int key, double *value) {
    size_t index = hash_function(map, key);
    for (size_t i = 0; i < map->max_buckets; i++) {
        size_t probe_index = (index + i) % map->max_buckets;
        if (map->buckets[probe_index].is_occupied && map->buckets[probe_index].key == key) {
            *value = map->buckets[probe_index].value;
            return true;
        }
        if (!map->buckets[probe_index].is_occupied) {
            break; // Key not found
        }
    }
    return false;
}
