#include "string_stream.h"

#include "kernel_prototypes.h"

#include <stdarg.h>
#include <stdio.h>
#include <string.h>

void string_stream_init(StringStream *s)
{
    s->length = 0; /* Does not include the terminating null character */
    s->capacity = 65536;
    s->buffer = (char*)my_malloc(s->capacity);
    s->buffer[0] = '\0';
}

void string_stream_sprintf(StringStream *s, const char * format, ...)
{
    if (s->buffer == NULL)
    {
        return;
    }

    va_list args;
    while (1)
    {
        const size_t capacity_left = s->capacity - s->length;
        va_start(args, format);
        const int needed = vsnprintf(
            s->buffer + s->length, capacity_left,
            format, args);
        va_end(args);

        if (needed < 0)
        {
            /* Formatting error */
            my_free(s->buffer);
            s->buffer = NULL;
            break;
        }

        /* Account for terminating null character. */
        if ((size_t)needed + 1 <= capacity_left)
        {
            s->length += needed;
            break;
        }

        s->capacity *= 2;
        char * const new_buffer = (char *)my_malloc(s->capacity);
        if (new_buffer)
        {
            memcpy(new_buffer, s->buffer, s->length + 1);
        }
        my_free(s->buffer);
        s->buffer = new_buffer;
    }
}
