#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <sys/resource.h>

/* Count allocation calls and requested bytes emitted by the Fortran
 * executable and its static library between allocation_begin and
 * allocation_end. Shared Fortran runtime internals are outside this
 * measurement. Requested bytes of realloc count in full; freed bytes
 * are not subtracted, so the byte total is a request total, not a
 * high-water mark. The high-water mark of the process is VmHWM of its
 * own address space: getrusage keeps the mark of the image before exec,
 * which floors small programs at the launching interpreter's size. */
void *__real_malloc(size_t size);
void *__real_calloc(size_t count, size_t size);
void *__real_realloc(void *pointer, size_t size);

static uint64_t allocation_calls;
static uint64_t allocation_bytes;
static int count_allocations;

void allocation_begin(void)
{
    allocation_calls = 0;
    allocation_bytes = 0;
    count_allocations = 1;
}

void allocation_end(int64_t *calls, int64_t *bytes)
{
    count_allocations = 0;
    *calls = (int64_t)allocation_calls;
    *bytes = (int64_t)allocation_bytes;
}

int64_t peak_rss_kilobytes(void)
{
    char line[256];
    long kilobytes = -1;
    FILE *status = fopen("/proc/self/status", "r");
    if (status) {
        while (fgets(line, sizeof line, status)) {
            if (strncmp(line, "VmHWM:", 6) == 0) {
                sscanf(line + 6, "%ld", &kilobytes);
                break;
            }
        }
        fclose(status);
    }
    if (kilobytes < 0) {
        struct rusage usage;
        getrusage(RUSAGE_SELF, &usage);
        kilobytes = usage.ru_maxrss;
    }
    return (int64_t)kilobytes;
}

void *__wrap_malloc(size_t size)
{
    if (count_allocations) {
        ++allocation_calls;
        allocation_bytes += size;
    }
    return __real_malloc(size);
}

void *__wrap_calloc(size_t count, size_t size)
{
    if (count_allocations) {
        ++allocation_calls;
        allocation_bytes += count * size;
    }
    return __real_calloc(count, size);
}

void *__wrap_realloc(void *pointer, size_t size)
{
    if (count_allocations) {
        ++allocation_calls;
        allocation_bytes += size;
    }
    return __real_realloc(pointer, size);
}
