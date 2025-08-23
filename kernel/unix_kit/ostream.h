/*
 * ostream.h
 *
 * An abstract stream that we can write to using printf (similar to
 * std::ostream in C++).
 *
 * It can be baked by a file or buffer that is automatically re-allocated
 * as needed (somewhat inspired by std::ofstream and std::stringstream).
 *
 * Example:
 *
 * void MyWrite(OStream * stream)
 * {
 *     for (int i = 0; i < 10000; i++ )
 *     {
 *         ostream_printf(stream, "%d\n", i);
 *     }
 * }
 *
 * void MyWriteToFile(FILE * fp)
 * {
 *     OStream stream;
 *     ofstream_init(&stream, fp);
 *     MyWrite(stream);
 * }
 *
 * char * MyString()
 * {
 *     OStream stream;
 *     string_stream_init(&stream);
 *     MyWrite(stream);
 *     return stream.buffer;
 * }
 *
 * Note: Ownership is transferred to callee of MyString. It has to call
 * my_free().
 *
 */

#ifndef _ostream_
#define _ostream_

#include <stdio.h>

typedef struct OStream OStream;

void string_stream_init(OStream *stream);
void ofstream_init(OStream *stream, FILE * file);
void ostream_printf(OStream *stream, const char * format, ...);

/************************************************************************/

/*
 * Implementation detail
 */
struct OStream
{
    FILE * file;

    char * buffer;
    size_t length;
    size_t capacity;
};

#endif
