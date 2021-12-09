/*
 *  words.c
 *
 *  This file provides functions to perform basic operations on
 *  words in a group. Such a word is given as a null-terminated array of
 *  integers. A positive integer g stands for the g-th generator and
 *  a negative integer -g for the inverse of the g-th generator.
 *
 *  The file provides the function
 *
 *      int *invert_group_word(int *word);
 *
 *  which reverses the order and the flips the signs, as well as the
 *  function
 *
 *      int *concat_group_words(int *word0, int *word1);
 *
 *  which concatenates two words and reduces the result by
 *  cancelling inverses that appear at the end of the first
 *  and beginning of the second word as much as possible.
 *
 *  The file also provides
 *
 *      int *copy_group_word(int *word);
 *
 *  to just make a copy of a null-terminated array of integers.
 */

/*
 * File added by MG 2022-01-24
 */

#include "kernel.h"
#include "kernel_namespace.h"

static int length_of_group_word(int *word);

/*
 * Compute length of null-terminated word, not including null terminator,
 * similar to strlen.
 */
static int length_of_group_word(
    int *word)
{
    int l = 0;
    while (*(word + l) != 0) l++;
    return l;
}

int* concat_group_words(
    int *word0,
    int *word1)
{
    int *prod;
    int *pos;
    int *word0_last;

    if (word0 == NULL)
        return NULL;
    if (word1 == NULL)
        return NULL;

    word0_last = word0 + length_of_group_word(word0) - 1;
    while (word0 <= word0_last && *word0_last == -*word1) {
        word0_last--; word1++;
    }

    prod = NEW_ARRAY(word0_last - word0 + length_of_group_word(word1) + 2, int);
    pos = prod;
    
    while (word0 <= word0_last)
        *(pos++) = *(word0++);

    while (*word1 != 0)
        *(pos++) = *(word1++);

    *pos = 0;
    
    return prod;
}

int* invert_group_word(
    int * word)
{
    int l;
    int *inv;
    int *back;

    if (word == NULL)
        return NULL;

    l = length_of_group_word(word);
    inv = NEW_ARRAY(l + 1, int);
    back = inv + l;
    *(back--) = 0;

    while (*word != 0)
        *(back--) = -*(word++);

    return inv;
}

int * copy_group_word(
    int *word)
{
    int *new_word;
    int *pos;

    if (word == NULL) {
        return NULL;
    }

    new_word = NEW_ARRAY(length_of_group_word(word) + 1, int);
    pos = new_word;
    while (*word != 0)
        *(pos++) = *(word++);

    *pos = 0;

    return new_word;
}
 
#include "end_namespace.h"
