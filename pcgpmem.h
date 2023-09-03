#ifndef PCGPMEM_H
#define PCGPMEM_H

#include <cstdint>
#include <cstddef>

#define PCGPMEM_ALIGN (2 * sizeof(void *))

typedef struct {
	void *start;
	void *end;
	void *ptr;
} Arena;

void arenaInit(Arena *a, void *mem, size_t size);
void *arenaAlloc(Arena *a, size_t size);
void arenaFree(Arena *a);

#endif
