#include "parse_util.h"

#include <stdlib.h>
#include <string.h>

char *parse_line(char **text)
{
    char * result;
    char c;

    if (*text == NULL) {
        return 0;
    }
    
    result = *text;

    while (**text && (**text != '\r') && (**text != '\n')) {
        (*text) ++;
    }

    c = **text;

    if (!c) {
        *text = NULL;
	return result;
    }

    **text = 0;
    (*text)++;

    if ( (c == '\r' && **text == '\n') ||
 	 (c == '\n' && **text == '\r')) {
        (*text)++;
    }
    
    return result;
}

static int is_line_empty(char *line)
{
    while(line) {
        if (*line != ' ' && *line != '\t') {
            return 0;
        }
        line++;
    }
    return 1;
}

char *parse_line_skipping_empty_lines(char **text)
{
    char * line;
    while ((line = parse_line(text)) && is_line_empty(line)) {
    }

    return line;
}

char *parse_token(char **line)
{
    char * result;

    if (*line == NULL) {
        return 0;
    }
    
    while (**line && (**line == ' ' || **line == '\t')) {
      (*line)++;
    }

    if (!**line) {
        return 0;
    }
    
    result = *line;

    while (**line && **line != ' ' && **line != '\t') {
      (*line)++;
    }
    
    if (**line) {
        **line = 0;
	(*line)++;
    }
    
    return result;
}

char *parse_token_next_non_empty_line(char **text, char **line)
{
    char * token;
    do {
        *line = parse_line(text);
        token = NULL;
    } while (*line && !(token = parse_token(line)));

    return token;
}

char *my_strdup(char *s)
{
    char *result = (char *) malloc(strlen(s) + 1);
    if (result != NULL) {
        strcpy(result, s);
    }
    return result;
}
