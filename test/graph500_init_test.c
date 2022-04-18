#include <stdio.h>
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <cmocka.h>
#include <aml.h>

/* A test case that does nothing and succeeds. */
static void graph500_init(void **state) {
    printf("Welcome to graph 500\n");
}

int main()
{
    const struct CMUnitTest tests[] = {
        cmocka_unit_test(graph500_init),
    };

    return cmocka_run_group_tests(tests, NULL, NULL);
}