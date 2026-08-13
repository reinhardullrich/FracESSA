#include <gtest/gtest.h>

#include <fracessa/fraction.hpp>

#include <sstream>
#include <utility>

using namespace fracessa::numeric;

/* Exact public-output value semantics around FLINT-backed storage. */

TEST(FractionTest, DefaultAndExplicitZero) {
    fraction value;
    EXPECT_EQ(value.to_string(), "0");
    EXPECT_EQ(value.to_dbl(), 0.0);

    value = fraction(3, 4);
    value.set_zero();
    EXPECT_EQ(value.to_string(), "0");
}

TEST(FractionTest, NumericConstructionCanonicalizes) {
    EXPECT_EQ(fraction(4, 8), fraction(1, 2));
    EXPECT_EQ(fraction(1, -2).to_string(), "-1/2");
    EXPECT_EQ(fraction(-1, -2).to_string(), "1/2");
    EXPECT_EQ(fraction(100LL, 4LL).to_string(), "25");
}

TEST(FractionTest, CopiesAndMovesValues) {
    const fraction original(3, 7);
    fraction copy(original);
    EXPECT_EQ(copy, original);

    fraction moved(std::move(copy));
    EXPECT_EQ(moved, original);

    fraction assigned;
    assigned = original;
    EXPECT_EQ(assigned, original);

    fraction move_assigned;
    move_assigned = std::move(assigned);
    EXPECT_EQ(move_assigned, original);
}

TEST(FractionTest, MaterializesIntegerRatio) {
    const integer numerator(-6);
    const integer denominator(8);
    fraction value;
    value.set_ratio(numerator, denominator);

    EXPECT_EQ(value.to_string(), "-3/4");
    EXPECT_DOUBLE_EQ(value.to_dbl(), -0.75);
}

TEST(FractionTest, StreamsCanonicalText) {
    std::ostringstream output;
    output << fraction(5, 6);
    EXPECT_EQ(output.str(), "5/6");
}
