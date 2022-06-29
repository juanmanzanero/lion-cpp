#include "lion/foundation/lion_exception.h"
#include "gtest/gtest.h"

TEST(lion_exception_test, catch_as_lion_exception)
{
    try
    {
        throw lion_exception("exception test");
        FAIL();
    }
    catch(lion_exception& ex)
    {
        EXPECT_EQ(std::string(ex.what()), std::string("exception test"));
        SUCCEED();
    }
    catch(...)
    {
        FAIL();
    }
}


TEST(lion_exception_test, catch_as_runtime_error)
{
    try
    {
        throw lion_exception("exception test");
        FAIL();
    }
    catch(lion_exception& ex)
    {
        EXPECT_EQ(std::string(ex.what()), std::string("exception test"));
        SUCCEED();
    }
    catch(...)
    {
        FAIL();
    }
}

TEST(lion_exception_test, catch_as_exception)
{
    try
    {
        throw lion_exception("exception test");
        FAIL();
    }
    catch(std::exception& ex)
    {
        EXPECT_EQ(std::string(ex.what()), std::string("exception test"));
        SUCCEED();
    }
    catch(...)
    {
        FAIL();
    }
}
