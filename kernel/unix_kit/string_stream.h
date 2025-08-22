/*
 *  string_stream.h
 *
 * Provides a safe version of sprintf, automatically re-allocating the buffer
 * if there is insufficient space.
 *
 * Example:
 *
 * StringStream ss;
 * string_stream_init(&ss);
 * for (int i = 0; i < 10000; i++)
 * {
 *     string_stream_sprintf(&ss, "%d\n", i);
 * }
 * printf("%s\n", ss.buffer);
 * my_free(ss.buffer);
 *
 */

#ifndef _string_stream_
#define _string_stream_

#include <stddef.h>

typedef struct StringStream StringStream;

void string_stream_init(StringStream *s);
void string_stream_sprintf(StringStream *s, const char * format, ...);

/************************************************************************/

/*
 * Implementation detail
 */
struct StringStream
{
    char * buffer;
    size_t length;
    size_t capacity;
};

#endif
