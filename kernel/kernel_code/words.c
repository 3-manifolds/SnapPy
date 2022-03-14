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

    /*
     * Look at last letter of first word.
     */
    word0_last = word0 + length_of_group_word(word0) - 1;

    /*
     * Eat last letter of first word and first letter of second word
     * while one is the inverse (the negative) of the other.
     */
    while (word0 <= word0_last && *word0_last == -*word1) {
        word0_last--; word1++;
    }

    /*
     * Allocate memory for result.
     *
     * Note that word0_last points to the byte before word0 if word0
     * is empty (or was completely cancelled).
     */
    prod = NEW_ARRAY((word0_last - word0 + 1) + length_of_group_word(word1) + 1,
                     int);

    /*
     * Copy each word into result.
     */
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

    /*
     * Allocate memory for result.
     */
    l = length_of_group_word(word);
    inv = NEW_ARRAY(l + 1, int);

    /*
     * Fill the result backwards starting with null-terminator.
     */
    back = inv + l;
    *(back--) = 0;

    /*
     * Copy into result backward.
     */
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

    /*
     * Allocate memory for result.
     */
    new_word = NEW_ARRAY(length_of_group_word(word) + 1, int);

    /*
     * And copy into it.
     */
    pos = new_word;
    while (*word != 0)
        *(pos++) = *(word++);

    *pos = 0;

    return new_word;
}
 
#include "end_namespace.h"
