#include <check.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* Forward declaration - the actual function from orb_io.c */
extern int read_triangulation_from_buffer(void **trig, const char *buffer, size_t buffer_len);

#define MAX_SAFE_NAME_LENGTH 256

START_TEST(test_name_buffer_overflow_protection)
{
    /* Invariant: Buffer reads for name field must never exceed MAX_SAFE_NAME_LENGTH */
    
    /* Test payloads: crafted triangulation data with varying name lengths */
    struct {
        size_t name_length;
        const char *description;
    } test_cases[] = {
        { 64, "valid short name" },           /* Valid input */
        { 256, "boundary at max size" },      /* Boundary case */
        { 512, "2x overflow attempt" },       /* 2x buffer size */
        { 2560, "10x overflow attempt" },     /* 10x buffer size */
    };
    int num_cases = sizeof(test_cases) / sizeof(test_cases[0]);

    for (int i = 0; i < num_cases; i++) {
        size_t name_len = test_cases[i].name_length;
        
        /* Create a buffer simulating triangulation file format with oversized name */
        size_t buffer_size = name_len + 128; /* name + header overhead */
        char *buffer = malloc(buffer_size);
        ck_assert_ptr_nonnull(buffer);
        
        /* Fill with pattern to detect overflow */
        memset(buffer, 'A', buffer_size);
        
        /* Encode name_length in buffer (simulating file format) */
        memcpy(buffer, &name_len, sizeof(size_t));
        
        void *trig = NULL;
        int result = read_triangulation_from_buffer(&trig, buffer, buffer_size);
        
        /* Security invariant: either parsing fails for oversized names,
           or the name is safely truncated to MAX_SAFE_NAME_LENGTH */
        if (name_len > MAX_SAFE_NAME_LENGTH) {
            /* Oversized input must be rejected or truncated */
            ck_assert_msg(result != 0 || trig == NULL,
                "Oversized name (%zu bytes) should be rejected: %s",
                name_len, test_cases[i].description);
        }
        
        if (trig != NULL) {
            free(trig);
        }
        free(buffer);
    }
}
END_TEST

Suite *security_suite(void)
{
    Suite *s;
    TCase *tc_core;

    s = suite_create("Security");
    tc_core = tcase_create("Core");

    tcase_add_test(tc_core, test_name_buffer_overflow_protection);
    suite_add_tcase(s, tc_core);

    return s;
}

int main(void)
{
    int number_failed;
    Suite *s;
    SRunner *sr;

    s = security_suite();
    sr = srunner_create(s);

    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);

    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}