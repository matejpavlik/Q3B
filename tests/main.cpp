#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "../lib/Solver.cpp"

void InitSylvan()
{
    using namespace sylvan;
    lace_init(1, 0);
    lace_startup(0, NULL, NULL);
    sylvan_set_sizes(1LL<<23, 1LL<<26, 1LL<<23, 1LL<<26);
    sylvan_init_package();
    sylvan_init_bdd();
}

Result SolveWithoutApprox(std::string filename)
{
    InitSylvan();
    Config config;
    config.propagateUnconstrained = true;
    config.approximationMethod = VARIABLES;
    Solver solver(config);

    z3::context ctx;
    z3::solver s(ctx);
    s.from_file(filename.c_str());
    z3::expr expr = mk_and(s.assertions());

    return solver.Solve(expr);
}

Result SolveWithVariableApprox(std::string filename, Approximation approx = NO_APPROXIMATION)
{
    InitSylvan();
    Config config;
    config.propagateUnconstrained = true;
    config.approximationMethod = VARIABLES;
    Solver solver(config);

    z3::context ctx;
    z3::solver s(ctx);
    s.from_file(filename.c_str());
    z3::expr expr = mk_and(s.assertions());

    return solver.Solve(expr, approx);
}

Result SolveWithOperationsLimitApprox(std::string filename, Approximation approx = NO_APPROXIMATION, int precision = 0)
{
    InitSylvan();
    Config config;
    config.propagateUnconstrained = true;
    config.approximationMethod = OPERATIONS;
    config.checkModels = true;
    Solver solver(config);

    z3::context ctx;
    z3::solver s(ctx);
    s.from_file(filename.c_str());
    z3::expr expr = mk_and(s.assertions());

    return solver.Solve(expr, approx, precision);
}

Result SolveWithBothLimitApprox(std::string filename, Approximation approx = NO_APPROXIMATION, int precision = 0)
{
    InitSylvan();
    Config config;
    config.propagateUnconstrained = true;
    config.approximationMethod = BOTH;
    config.checkModels = true;
    Solver solver(config);

    z3::context ctx;
    z3::solver s(ctx);
    s.from_file(filename.c_str());
    z3::expr expr = mk_and(s.assertions());

    return solver.Solve(expr, approx, precision);
}

TEST_CASE( "Without approximations", "[noapprox]" )
{
    REQUIRE( SolveWithoutApprox("../tests/data/AR-fixpoint-1.smt2") == UNSAT );
    REQUIRE( SolveWithoutApprox("../tests/data/cache-coherence-2-fixpoint-1.smt2") == UNSAT );
    REQUIRE( SolveWithoutApprox("../tests/data/itc-b13-fixpoint-3.smt2") == SAT );
    REQUIRE( SolveWithoutApprox("../tests/data/Fibonacci01_true-unreach-call_true-no-overflow.c_905.smt2") == SAT );
    REQUIRE( SolveWithoutApprox("../tests/data/nlzbe008.smt2") == UNSAT );
    REQUIRE( SolveWithoutApprox("../tests/data/falseAndFalse.smt2") == UNSAT );
    REQUIRE( SolveWithoutApprox("../tests/data/bvshl0.smt2") == SAT );
    REQUIRE( SolveWithoutApprox("../tests/data/check_bvsge_bvashr0_16bit.smt2") == UNSAT );
    REQUIRE( SolveWithoutApprox("../tests/data/magnetic_field-node118398.smt2") == UNSAT );
    REQUIRE( SolveWithoutApprox("../tests/data/unconstrainedMulVar.smt2") == UNSAT );
    REQUIRE( SolveWithoutApprox("../tests/data/check_bvsle_bvmul_8bit.smt2") == UNSAT );
    REQUIRE( SolveWithoutApprox("../tests/data/unconstrainedMulConst.smt2") == SAT );
    REQUIRE( SolveWithoutApprox("../tests/data/check_eq_bvconcat0_2_64bit.smt2") == UNSAT );
    REQUIRE( SolveWithoutApprox("../tests/data/002.smt2") == UNSAT );
}

TEST_CASE( "With variable approximations", "[variableapprox]" )
{
    REQUIRE( SolveWithVariableApprox("../tests/data/audio_ac97_wavepcistream2.cpp.smt2", OVERAPPROXIMATION) == UNSAT );
    REQUIRE( SolveWithVariableApprox("../tests/data/jain_7_true-unreach-call_true-no-overflow.i_61.smt2", OVERAPPROXIMATION) == UNSAT );
    REQUIRE( SolveWithVariableApprox("../tests/data/RNDPRE_3_48.smt2", UNDERAPPROXIMATION) == SAT );
    REQUIRE( SolveWithVariableApprox("../tests/data/ETCS-essentials-node3023.smt2", UNDERAPPROXIMATION) == SAT );
    REQUIRE( SolveWithVariableApprox("../tests/data/accelerating-node2100.smt2", UNDERAPPROXIMATION) == SAT );
    REQUIRE( SolveWithVariableApprox("../tests/data/binary_driver-2007-10-09-node11383.smt2", UNDERAPPROXIMATION) == SAT );
    REQUIRE( SolveWithVariableApprox("../tests/data/ARI118=1.smt2", OVERAPPROXIMATION) == UNSAT );
}

TEST_CASE( "With bothLimit approximations", "[bothlimitapprox]" )
{
    REQUIRE( SolveWithBothLimitApprox("../tests/data/RNDPRE_4_42.smt2", OVERAPPROXIMATION) == UNSAT );
    REQUIRE( SolveWithBothLimitApprox("../tests/data/RND_6_4.smt2", UNDERAPPROXIMATION) == SAT );
    REQUIRE( SolveWithBothLimitApprox("../tests/data/jain_7_true-unreach-call_true-no-overflow.i_61.smt2", OVERAPPROXIMATION) == UNSAT );

    //correct model returned by an overapproximation
    REQUIRE( SolveWithBothLimitApprox("../tests/data/007.smt2", OVERAPPROXIMATION) == SAT );
    REQUIRE( SolveWithBothLimitApprox("../tests/data/sum02_true-unreach-call_true-no-overflow.i_375.smt2", OVERAPPROXIMATION) == SAT );
    REQUIRE( SolveWithBothLimitApprox("../tests/data/check_bvsgt_bvudiv1_8bit.smt2", UNDERAPPROXIMATION) != SAT );
    REQUIRE( SolveWithBothLimitApprox("../tests/data/bvurem_approx.smt2", UNDERAPPROXIMATION, 1) != SAT );
}

TEST_CASE( "With bothLimit approximations -- term introducer ", "[bothlimitapprox-ti]" )
{
    REQUIRE( SolveWithBothLimitApprox("../tests/data/intersection-example-onelane.proof-node1469.smt2", OVERAPPROXIMATION) == UNSAT );
}

TEST_CASE( "With operation approximations -- ite ", "[opapproxlimit-ite]" )
{
    REQUIRE( SolveWithOperationsLimitApprox("../tests/data/iteApprox.smt2", UNDERAPPROXIMATION, 1) != UNSAT );
}

TEST_CASE( "SMT-COMP 2018", "[smtcomp18]" )
{
    REQUIRE( SolveWithVariableApprox( "../tests/data/smtcomp18/01.smt2", UNDERAPPROXIMATION ) != SAT );
    REQUIRE( SolveWithBothLimitApprox( "../tests/data/smtcomp18/02.smt2", OVERAPPROXIMATION, 1 ) != UNSAT );
}
