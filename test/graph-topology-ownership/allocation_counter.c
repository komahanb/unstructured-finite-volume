#include <stddef.h>
#include <stdint.h>

/* Count calls emitted by the Fortran executable and its static library.
 * Shared Fortran runtime internals are outside this measurement. */
void *__real_malloc(size_t size);
void *__real_calloc(size_t count, size_t size);
void *__real_realloc(void *pointer, size_t size);

static uint64_t allocation_count;
static int count_allocations;

void allocation_begin(void)
{
    allocation_count = 0;
    count_allocations = 1;
}

int64_t allocation_end(void)
{
    count_allocations = 0;
    return (int64_t)allocation_count;
}

void *__wrap_malloc(size_t size)
{
    if (count_allocations) ++allocation_count;
    return __real_malloc(size);
}

void *__wrap_calloc(size_t count, size_t size)
{
    if (count_allocations) ++allocation_count;
    return __real_calloc(count, size);
}

void *__wrap_realloc(void *pointer, size_t size)
{
    if (count_allocations) ++allocation_count;
    return __real_realloc(pointer, size);
}
