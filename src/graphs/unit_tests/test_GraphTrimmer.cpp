//
// Created by sylwester on 3/18/21.
//

#include "graphs/GraphTrimmer.h"
#include "gtest/gtest.h"


// TESTS
class GraphTrimmerFixture : public ::testing::Test {
protected:
    virtual void TearDown() {
        clog << "TearDown, V.clear()" << endl;
    }

    virtual void SetUp() {
        clog << "SetUp V = ..." << endl;
        V = {
                {3,7}, // 0
                {5,8}, // 1
                {18,19}, // 2
                {0,4}, // 3
                {3,7}, // 4
                {1,8}, // 5
                {19,14}, // 6
                {0,4}, // 7
                {1,5}, // 8
                {15,18}, // 9
                {14,13}, // 10
                {13,12}, // 11
                {11}, // 12
                {10,11}, // 13
                {18,10,6}, // 14
                {16,17,9}, // 15
                {15}, // 16
                {15}, // 17
                {9,2,14}, // 18
                {2,6} // 19
        };
    }

    VVI V;
};


TEST_F(GraphTrimmerFixture, findMaximalPathsAndCycles) {
    GraphTrimmer trm(V);
    auto pth_cyc = trm.findMaximalPathsAndCycles();

    EXPECT_EQ( pth_cyc.size(), 7 );
    EXPECT_EQ( pth_cyc[0].size(), 4 );
    EXPECT_EQ( pth_cyc[1].size(), 3 );
    EXPECT_EQ( pth_cyc[2].size(), 3 );
    EXPECT_EQ( pth_cyc[3].size(), 1 );
    EXPECT_EQ( pth_cyc[4].size(), 4 );
    EXPECT_EQ( pth_cyc[5].size(), 1 );
    EXPECT_EQ( pth_cyc[6].size(), 1 );
}

TEST_F(GraphTrimmerFixture, findMaximalPathsAndCycles2) {
    GraphTrimmer trm(V);
    auto pth_cyc = trm.findMaximalPathsAndCycles();

    EXPECT_EQ( pth_cyc.size(), 7 );
    EXPECT_EQ( pth_cyc[0].size(), 4 );
    EXPECT_EQ( pth_cyc[1].size(), 3 );
    EXPECT_EQ( pth_cyc[2].size(), 3 );
    EXPECT_EQ( pth_cyc[3].size(), 1 );
    EXPECT_EQ( pth_cyc[4].size(), 4 );
    EXPECT_EQ( pth_cyc[5].size(), 1 );
    EXPECT_EQ( pth_cyc[6].size(), 1 );
}