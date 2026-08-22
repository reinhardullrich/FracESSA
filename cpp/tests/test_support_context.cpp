#include <gtest/gtest.h>

#include <vector>

#include <fracessa/support.hpp>

using fracessa::support::Support;
using fracessa::support::SupportContext;

TEST(SupportContextTest, SmallAndLargeRepresentationsHaveTheSameSemantics)
{
    for (const size_t dimension : {20u, 64u, 65u, 129u}) {
        SupportContext context(dimension);
        Support support = context.make();
        const std::vector<size_t> expected_indices{0, dimension / 2, dimension - 1};
        for (const size_t position : expected_indices) context.set(support, position);

        EXPECT_EQ(context.count(support), expected_indices.size());
        std::vector<size_t> indices;
        context.extract_set_indices(support, indices);
        EXPECT_EQ(indices, expected_indices);

        Support clone = context.clone(support);
        context.reset(clone, dimension / 2);
        EXPECT_TRUE(context.is_subset_of(clone, support));
        EXPECT_FALSE(context.equal(clone, support));

        Support rotated = context.clone(support);
        for (size_t shift = 0; shift < dimension; ++shift) context.rotate_one_right(rotated);
        EXPECT_TRUE(context.equal(rotated, support));

        Support reflected = context.clone(support);
        context.reflect(reflected);
        context.reflect(reflected);
        EXPECT_TRUE(context.equal(reflected, support));

        Support complement = context.make();
        context.set_all(complement);
        context.subtract(complement, support);
        EXPECT_EQ(context.count(complement), dimension - expected_indices.size());
    }
}

TEST(SupportContextTest, RejectsZeroDimension)
{
    EXPECT_THROW(SupportContext(0), std::invalid_argument);
}
