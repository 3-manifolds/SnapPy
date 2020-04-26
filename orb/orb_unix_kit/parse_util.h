#ifndef _parsing_util_
#define _parsing_util_

extern char *parse_line(char **text);

extern char *parse_line_skipping_empty_lines(char **text);

extern char *parse_token(char **line);

extern char *parse_token_next_non_empty_line(char **text, char **line);

extern char *my_strdup(char *s);

#endif
