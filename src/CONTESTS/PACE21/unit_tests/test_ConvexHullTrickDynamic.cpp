//
// Created by sylwester on 4/1/21.
//

#include <utils/RandomNumberGenerators.h>
#include <Constants.h>
#include <utils/TimeMeasurer.h>
#include "CONTESTS/PACE21/heur/ConvexHullTrickDynamic.h"
#include "gtest/gtest.h"

TEST( ConvexHullTrick, test_maximize_1 ){

    bool maximize = true;
    ConvexHullTrickDynamic hull( maximize );
    hull.insert_line( 3,-3,0 );
    hull.insert_line( 0,0,1 );
    hull.insert_line( -1,3,2 );

    ASSERT_EQ( hull.get_line_id( -100 ), 2  );
    ASSERT_EQ( hull.get_line_id( -1 ), 2  );
    ASSERT_EQ( hull.get_line_id( 0 ), 2  );

    ASSERT_TRUE( hull.get_line_id( 1 ) == 2 );

    ASSERT_EQ( hull.get_line_id( 2 ), 0  );

    ASSERT_TRUE( hull.get_line_id( 3 ) == 0 );

    ASSERT_EQ( hull.get_line_id( 4 ), 0  );
    ASSERT_EQ( hull.get_line_id( 5 ), 0  );
    ASSERT_EQ( hull.get_line_id( 100 ), 0  );
}

TEST( ConvexHullTrick, test_minimize_1 ){

    bool maximize = false;
    ConvexHullTrickDynamic hull( maximize );
    hull.insert_line( 3,-3,0 );
    hull.insert_line( 0,0,1 );
    hull.insert_line( -1,3,2 );

    ASSERT_EQ( hull.get_line_id( -100 ), 0  );
    ASSERT_EQ( hull.get_line_id( -1 ), 0  );
    ASSERT_EQ( hull.get_line_id( 0 ), 0  );

    ASSERT_TRUE( hull.get_line_id( 1 ) == 0 || hull.get_line_id( 1 ) == 1 );

    ASSERT_EQ( hull.get_line_id( 2 ), 1  );

    ASSERT_TRUE( hull.get_line_id( 3 ) == 1 || hull.get_line_id( 3 ) == 2 );

    ASSERT_EQ( hull.get_line_id( 4 ), 2  );
    ASSERT_EQ( hull.get_line_id( 5 ), 2  );
    ASSERT_EQ( hull.get_line_id( 100 ), 2  );
}

TEST( ConvexHullTrick, test_large ){

    // there will be roughly 2/N lines added
    // there will be roughly 2*R queries

//    double R = 5'000;
//    double N = 0.0003; // large test - large differences in performance

    double R = 200;
    double N = 0.005; // small test - almost no difference in performance

    cout << fixed;
    cout.precision(2);

    int REPS = 100;
    while(REPS--) {
        for (int maxim = 0; maxim <= 1; maxim++) {

            bool maximize = maxim;
            ConvexHullTrickDynamic hull(maximize);

            vector<pair<double, double>> lines;

            for (double i = -R; i <= R; i += (N * R)) {
                double x = i;
                pair<double, double> line;
                if (x == 0) {
                    line = {0, R};
                    lines.push_back(line);
//                cout << "Line for i: " << i << "  |  " << lines.back().first << " * x  +  " << lines.back().second << endl;
                    continue;
                }

                double y = sqrt(1ll * R * R - x * x);
                if (y < 0.01) continue;

                double ctga = x / y;

                {
                    // mx + b = y
                    double m = -ctga; // tg(90+alpha) = -ctg(alpha)
                    double b = y - m * x;
                    line = {m, b};
                    lines.push_back(line);
                }

//            cout << "Line for i: " << i << "  |  " << lines.back().first << " * x  +  " << lines.back().second << endl;
            }

//        cout << "lines.size(): " << lines.size() << ", lines:" << endl << lines << endl;

            TimeMeasurer::startMeasurement("hull");
            for (int i = 0; i < lines.size(); i++) hull.insert_line(lines[i].first, lines[i].second, i);
            TimeMeasurer::stopMeasurement("hull");

            auto brute = [&](double x) {
                double m = 1e9;
                int id = -1;
                if (maximize) m = -1e9;
                for (int i = 0; i < lines.size(); i++) {
                    double val = lines[i].first * x + lines[i].second;
                    if (maximize && (val > m)) {
                        m = val;
                        id = i;
                    } else if (!maximize && (val < m)) {
                        m = val;
                        id = i;
                    }
                }
                return id;
            };

//        double STEP = N*R;
            double STEP = 1;
            unordered_set<int> zb;
            for (double i = -R; i <= R; i += STEP) {
//            cout << "i: " << i << endl;

                TimeMeasurer::startMeasurement("brute");
                auto brute_eval_id = brute(i);
                TimeMeasurer::stopMeasurement("brute");

                TimeMeasurer::startMeasurement("hull");
                auto hull_eval_id = hull.get_line_id(i);
                TimeMeasurer::stopMeasurement("hull");

                auto brute_val = lines[brute_eval_id].first * i + lines[brute_eval_id].second;
                auto hull_val =
                        lines[hull.get_line_id(hull_eval_id)].first * i + lines[hull.get_line_id(hull_eval_id)].second;

                if (brute_eval_id != hull_eval_id) {
//                    cout << "brute_eval_id: " << brute_eval_id << endl;
//                    cout << "hull_eval_id: " << hull_eval_id << endl;
//                    cout << "brute val: " << lines[brute_eval_id].first * i + lines[brute_eval_id].second << endl;
//                    cout << "hull val: " << lines[hull.get_line_id(hull_eval_id)].first * i +
//                                            lines[hull.get_line_id(hull_eval_id)].second << endl;
                    ASSERT_NEAR(brute_val, hull_val, 0.00000001);
                } else
                    ASSERT_EQ(brute_eval_id, hull_eval_id);

                zb.insert(hull.get_line_id(i));
                zb.insert(brute_eval_id);
            }

//            cout << "There are " << zb.size() << " different ids in zb" << endl;

//            cout << endl << endl;
        }
    }

    TimeMeasurer::writeAllMeasurements();
}