#include <gtsam/geometry/Pose2.h>
#include <gtsam/inference/Key.h>
#include <gtsam/slam/BetweenFactor.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/GaussNewtonOptimizer.h>
#include <gtsam/nonlinear/Marginals.h>
#include <gtsam/nonlinear/Values.h>

#include <iostream>

using namespace std;
using namespace gtsam;

int main()
{
    std::cout << "Simple factor graph..." << std::endl;
    // // 1. Create a factor graph container and add factors to it
    NonlinearFactorGraph graph;

    // // 2a. Add a prior on the first pose, setting it to the origin
    // // A prior factor consists of a mean and a noise model (covariance matrix)
    auto priorNoise = noiseModel::Diagonal::Sigmas(Vector3(0.3, 0.3, 0.1));
    // graph.addPrior(1, Pose2(0, 0, 0), priorNoise);

    // // For simplicity, we will use the same noise model for odometry and loop closures
    // auto model = noiseModel::Diagonal::Sigmas(Vector3(0.2, 0.2, 0.1));

    // // 2b. Add odometry factors
    // // Create odometry (Between) factors between consecutive poses
    // graph.emplace_shared<BetweenFactor<Pose2> >(1, 2, Pose2(2, 0, 0), model);
    // graph.emplace_shared<BetweenFactor<Pose2> >(2, 3, Pose2(2, 0, M_PI_2), model);
    // graph.emplace_shared<BetweenFactor<Pose2> >(3, 4, Pose2(2, 0, M_PI_2), model);
    // graph.emplace_shared<BetweenFactor<Pose2> >(4, 5, Pose2(2, 0, M_PI_2), model);

    // // 2c. Add the loop closure constraint
    // // This factor encodes the fact that we have returned to the same pose. In real systems,
    // // these constraints may be identified in many ways, such as appearance-based techniques
    // // with camera images. We will use another Between Factor to enforce this constraint:
    // graph.emplace_shared<BetweenFactor<Pose2> >(5, 2, Pose2(2, 0, M_PI_2), model);
    // graph.print("\nFactor Graph:\n");  // print

    // // 3. Create the data structure to hold the initialEstimate estimate to the solution
    // // For illustrative purposes, these have been deliberately set to incorrect values
    // Values initialEstimate;
    // initialEstimate.insert(1, Pose2(0.5, 0.0, 0.2));
    // initialEstimate.insert(2, Pose2(2.3, 0.1, -0.2));
    // initialEstimate.insert(3, Pose2(4.1, 0.1, M_PI_2));
    // initialEstimate.insert(4, Pose2(4.0, 2.0, M_PI));
    // initialEstimate.insert(5, Pose2(2.1, 2.1, -M_PI_2));
    // initialEstimate.print("\nInitial Estimate:\n");  // print

    // // 4. Optimize the initial values using a Gauss-Newton nonlinear optimizer
    // // The optimizer accepts an optional set of configuration parameters,
    // // controlling things like convergence criteria, the type of linear
    // // system solver to use, and the amount of information displayed during
    // // optimization. We will set a few parameters as a demonstration.
    // GaussNewtonParams parameters;
    // // Stop iterating once the change in error between steps is less than this value
    // parameters.relativeErrorTol = 1e-5;
    // // Do not perform more than N iteration steps
    // parameters.maxIterations = 100;
    // // Create the optimizer ...
    // GaussNewtonOptimizer optimizer(graph, initialEstimate, parameters);
    // // ... and optimize
    // Values result = optimizer.optimize();
    // result.print("Final Result:\n");

    // // 5. Calculate and print marginal covariances for all variables
    // cout.precision(3);
    // Marginals marginals(graph, result);
    // cout << "x1 covariance:\n" << marginals.marginalCovariance(1) << endl;
    // cout << "x2 covariance:\n" << marginals.marginalCovariance(2) << endl;
    // cout << "x3 covariance:\n" << marginals.marginalCovariance(3) << endl;
    // cout << "x4 covariance:\n" << marginals.marginalCovariance(4) << endl;
    // cout << "x5 covariance:\n" << marginals.marginalCovariance(5) << endl;
}