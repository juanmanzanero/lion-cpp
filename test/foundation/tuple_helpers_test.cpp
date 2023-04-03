#include <array>

#include "gtest/gtest.h"

#include "lion/foundation/tuple_helpers.h"


//
// Defines tests for the compile-time loops
// in header "lion/foundation/tuple_helpers.h".
//


constexpr std::size_t array_size = 100u;


TEST(tuple_helpers_test, index_for_total_traverse_test)
{
    std::array<double, array_size> arr_compiletime{};
    auto arr_runtime = arr_compiletime;

    index_for<array_size>(
        [&arr_compiletime](auto i) { arr_compiletime[i] += 10. * i; });

    for (auto i = 0u; i < array_size; ++i) {
        arr_runtime[i] += 10. * i;
    }

    EXPECT_EQ(arr_runtime, arr_compiletime);
}


TEST(tuple_helpers_test, index_for_zero_traverse_test)
{
    std::array<double, array_size> arr_compiletime{};
    auto arr_runtime = arr_compiletime;

    index_for<0u>(
        [&arr_compiletime](auto i) { arr_compiletime[i] += 10. * i; });

    EXPECT_EQ(arr_runtime, arr_compiletime);
}


TEST(tuple_helpers_test, index_for_partial_traverse_test)
{
    std::array<double, array_size> arr_compiletime{};
    auto arr_runtime = arr_compiletime;

    index_for<array_size / 2u, array_size>(
        [&arr_compiletime](auto i) { arr_compiletime[i] += 10. * i; });

    for (auto i = array_size / 2u; i < array_size; ++i) {
        arr_runtime[i] += 10. * i;
    }

    EXPECT_EQ(arr_runtime, arr_compiletime);
}


TEST(tuple_helpers_test, index_for_empty_traverse_test)
{
    std::array<double, array_size> arr_compiletime{};
    auto arr_runtime = arr_compiletime;

    index_for<array_size + 1u, array_size>(
        [&arr_compiletime](auto i) { arr_compiletime[i] += 10. * i; });

    EXPECT_EQ(arr_runtime, arr_compiletime);
}


TEST(tuple_helpers_test, tuple_for_test)
{
    std::array<double, array_size> arr_compiletime{};
    auto arr_runtime = arr_compiletime;

    std::size_t count = 0u;
    tuple_for(arr_compiletime, [&count](auto &aci) mutable { aci += 10. * count++; });

    count = 0u;
    for (auto &ari : arr_runtime) {
        ari += 10. * count++;
    }

    EXPECT_EQ(arr_runtime, arr_compiletime);
}


TEST(tuple_helpers_test, tuple_index_for_test)
{
    std::array<double, array_size> arr_compiletime{};
    auto arr_runtime = arr_compiletime;

    tuple_index_for(arr_compiletime, [](auto &aci, auto i) { aci += 10. * i; });

    for (auto i = 0u; i < array_size; ++i) {
        arr_runtime[i] += 10. * i;
    }

    EXPECT_EQ(arr_runtime, arr_compiletime);
}
