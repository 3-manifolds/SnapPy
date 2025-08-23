#include "ostream.h"

#include "kernel.h"

#include <stdarg.h>
#include <stdio.h>
#include <string.h>

void string_stream_init(OStream *s)
{
    s->file = NULL;

    s->length = 0; /* Does not include the terminating null character */
    s->capacity = 65536;
    s->buffer = (char*)my_malloc(s->capacity);
    s->buffer[0] = '\0';
}

void ofstream_init(OStream *s, FILE * file)
{
    s->file = file;

    s->length = 0;
    s->capacity = 0;
    s->buffer = NULL;
}

void ostream_printf(OStream *s, const char * format, ...)
{
    va_list args;

    if (s->file != NULL)
    {
        va_start(args, format);
        vfprintf(s->file, format, args);
        va_end(args);
	return;
    }
    
    while (1)
    {
        if (s->buffer == NULL)
        {
            break;
        }

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
