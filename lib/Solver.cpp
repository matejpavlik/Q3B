#include "Solver.h"
#include "ExprSimplifier.h"
#include "TermConstIntroducer.h"
#include "Logger.h"

#include <thread>
#include <functional>
#include <sstream>

std::mutex Solver::m;
std::mutex Solver::m_z3context;

Result Solver::getResult(z3::expr expr, Approximation approximation, int effectiveBitWidth)
{
    if (expr.is_const())
    {
	if (expr.is_app() && expr.decl().decl_kind() == Z3_OP_TRUE)
        {
            return SAT;
        }
	else if (expr.is_app() && expr.decl().decl_kind() == Z3_OP_FALSE)
        {
            return UNSAT;
        }
    }

    ExprToBDDTransformer transformer(expr.ctx(), expr, config);

    if (approximation == OVERAPPROXIMATION || approximation == UNDERAPPROXIMATION)
    {
	if (effectiveBitWidth == 0)
	{
	    if (approximation == OVERAPPROXIMATION)
	    {
		return runWithOverApproximations(transformer);
	    }
	    else
	    {
		return runWithUnderApproximations(transformer);
	    }
	}

        if (approximation == OVERAPPROXIMATION)
        {
            return runOverApproximation(transformer, effectiveBitWidth, abs(effectiveBitWidth));
        }
        else
        {
            return runUnderApproximation(transformer, effectiveBitWidth, abs(effectiveBitWidth));
        }
    }

    auto returned = transformer.Proccess();
    return returned.isZero() ? UNSAT : SAT;
}

Result Solver::Solve(z3::expr expr, Approximation approximation, int effectiveBitWidth)
{
    Logger::Log("Solver", "Simplifying formula.", 1);
    m_z3context.lock();
    ExprSimplifier simplifier(expr.ctx(), config.propagateUnconstrained);
    expr = simplifier.Simplify(expr);
    if (config.approximationMethod == OPERATIONS || config.approximationMethod == BOTH)
    {
	expr = simplifier.StripToplevelExistentials(expr);
    }

    m_z3context.unlock();

    bool negated = false;
    if (config.flipUniversalQuantifier &&
	expr.is_quantifier() &&
	Z3_is_quantifier_forall(expr.ctx(), (Z3_ast)expr) &&
	simplifier.isSentence(expr))
    {
	Logger::Log("Solver", "Negating universal formula.", 1);
	negated = true;
	expr = simplifier.negate(expr);
	expr = simplifier.PushNegations(expr);
	expr = simplifier.StripToplevelExistentials(expr);
	if (approximation == OVERAPPROXIMATION) approximation = UNDERAPPROXIMATION;
	else if (approximation == UNDERAPPROXIMATION) approximation = OVERAPPROXIMATION;
    }

    if (approximation == OVERAPPROXIMATION)
    {
	Logger::Log("Solver", "Introducing mul constants.", 1);
	TermConstIntroducer tci(expr.ctx());
	expr = tci.FlattenMul(expr);
    }

    Logger::Log("Solver", "Starting solver.", 1);
    auto result = getResult(expr, approximation, effectiveBitWidth);
    if (negated)
    {
	Logger::Log("Solver", "Flipping result of the negated formula.", 1);
	if (result == SAT) result = UNSAT;
	else if (result == UNSAT) result = SAT;
    }
    return result;
}

Result Solver::solverThread(z3::expr expr, Approximation approximation, int effectiveBitWidth)
{
    m_z3context.lock();
    z3::context ctx;
    auto translated = z3::to_expr(ctx, Z3_translate(expr.ctx(), expr, ctx));
    m_z3context.unlock();

    auto res = getResult(translated, approximation, effectiveBitWidth);

    if (res == SAT || res == UNSAT)
    {
	std::stringstream ss;

	if (approximation == NO_APPROXIMATION)
	{
	    Logger::Log("Solver", "Decided by the base solver", 1);
	}

	std::unique_lock<std::mutex> lk(m);
	resultComputed = true;
	result = res;
	doneCV.notify_one();
    }

    return res;
}

Result Solver::SolveParallel(z3::expr expr)
{
    Logger::Log("Solver", "Simplifying formula.", 1);
    ExprSimplifier simplifier(expr.ctx(), config.propagateUnconstrained);
    expr = simplifier.Simplify(expr);
    if (config.approximationMethod == OPERATIONS || config.approximationMethod == BOTH)
    {
	expr = simplifier.StripToplevelExistentials(expr);
    }

    if (expr.is_const())
    {
	if (expr.is_app() && expr.decl().decl_kind() == Z3_OP_TRUE)
        {
	    Logger::Log("Solver", "Solved by simplifications.", 1);
            return SAT;
        }
	else if (expr.is_app() && expr.decl().decl_kind() == Z3_OP_FALSE)
        {
	    Logger::Log("Solver", "Solved by simplifications.", 1);
            return UNSAT;
        }
    }

    bool negated = false;
    if (config.flipUniversalQuantifier &&
	expr.is_quantifier() &&
	Z3_is_quantifier_forall(expr.ctx(), (Z3_ast)expr) &&
	simplifier.isSentence(expr))
    {
	Logger::Log("Solver", "Negating universal formula.", 1);
	negated = true;
	expr = simplifier.negate(expr);
	expr = simplifier.PushNegations(expr);
	expr = simplifier.StripToplevelExistentials(expr);
    }

    Logger::Log("Solver", "Introducing mul constants.", 1);
    TermConstIntroducer tci(expr.ctx());
    auto overExpr = tci.FlattenMul(expr);

    Logger::Log("Solver", "Starting solver threads.", 1);
    auto main = std::thread( [this,expr] { solverThread(expr); } );
    main.detach();
    auto under = std::thread( [this,expr] { solverThread(expr, UNDERAPPROXIMATION); } );
    under.detach();
    auto over = std::thread( [this,overExpr] { solverThread(overExpr, OVERAPPROXIMATION); } );
    over.detach();

    std::unique_lock<std::mutex> lk(m);
    while (!resultComputed)
    {
	doneCV.wait(lk);
    }

    if (negated)
    {
	Logger::Log("Solver", "Flipping result of the negated formula.", 1);
	if (result == SAT) result = UNSAT;
	else if (result == UNSAT) result = SAT;
    }
    return result;
}

Result Solver::runOverApproximation(ExprToBDDTransformer &transformer, int bitWidth, int precision)
{
    transformer.setApproximationType(SIGN_EXTEND);

    std::stringstream ss;
    ss << "Trying bit-width " << bitWidth << ", precision " << precision;
    Logger::Log("Overapproximating solver", ss.str(), 5);

    auto returned = transformer.ProcessOverapproximation(bitWidth, precision);

    auto result = returned.upper.isZero() ? UNSAT : SAT;
    if (result == UNSAT)
    {
	Logger::Log("Solver", "Decided by overapproximation", 1);
	std::stringstream rss;
	rss << "Final bit-width " << bitWidth << ", precision " << precision;
	Logger::Log("Overapproximating solver", rss.str(), 1);
	return result;
    }

    transformer.PrintNecessaryValues(returned.upper);

    if (config.checkModels)
    {
	auto model = transformer.GetModel(returned.lower.isZero() ? returned.upper : returned.lower);

	m_z3context.lock();
	auto substituted = substituteModel(transformer.expression, model).simplify();
	m_z3context.unlock();

	if (substituted.hash() != transformer.expression.hash())
	{
	    Logger::Log("Overapproximating solver", "Validating model", 5);

	    Config validatingConfig;
	    validatingConfig.propagateUnconstrained = true;
	    validatingConfig.approximationMethod = config.approximationMethod;

	    Solver validatingSolver(validatingConfig);

	    if (validatingSolver.Solve(substituted, UNDERAPPROXIMATION, 1) == SAT)
	    {
		Logger::Log("Solver", "Decided by overapproximation", 1);
		std::stringstream rss;
		rss << "Final bit-width " << bitWidth << ", precision " << precision;
		Logger::Log("Overapproximating solver", rss.str(), 1);
		return SAT;
	    }
	}
    }

    return UNKNOWN;
}

Result Solver::runUnderApproximation(ExprToBDDTransformer &transformer, int bitWidth, int precision)
{
    transformer.setApproximationType(ZERO_EXTEND);

    std::stringstream ss;
    ss << "Trying bit-width " << bitWidth << ", precision " << precision;
    Logger::Log("Underapproximating solver", ss.str(), 5);

    auto returned = transformer.ProcessUnderapproximation(bitWidth, precision);
    auto result = returned.lower.isZero() ? UNSAT : SAT;

    if (result == SAT)
    {
	Logger::Log("Solver", "Decided by underapproximation", 1);
	std::stringstream rss;
	rss << "Final bit-width " << bitWidth << ", precision " << precision;
	Logger::Log("Underapproximating solver", rss.str(), 1);
	return result;
    }

    return UNKNOWN;
}

Result Solver::runWithApproximations(ExprToBDDTransformer &transformer, Approximation approximation)
{
    assert (approximation != NO_APPROXIMATION);

    auto runFunction = [this, &approximation](auto &transformer, int bitWidth, unsigned int precision){
	return (approximation == UNDERAPPROXIMATION) ?
	   runUnderApproximation(transformer, bitWidth, precision) :
	   runOverApproximation(transformer, bitWidth, precision);
    };

    if (config.approximationMethod == BOTH)
    {
	unsigned int prec = 1;
	unsigned int lastBW = 1;
	while (prec != 0)
	{
	    if (prec == 4 && approximation == OVERAPPROXIMATION)
	    {
		Result approxResult = runFunction(transformer, 32, 2);
		if (approxResult != UNKNOWN)
		{
		    return approxResult;
		}
	    }

	    if (lastBW == 1)
	    {
		Result approxResult = runFunction(transformer, lastBW, prec);
		if (approxResult != UNKNOWN)
		{
		    return approxResult;
		}

		bool approxHappened = transformer.OperationApproximationHappened();

		if (approxHappened || transformer.OperationApproximationHappened())
		{
		    prec *= 4;
		    continue;
		}

		approxResult = runFunction(transformer, -1, prec);
		if (approxResult != UNKNOWN)
		{
		    return approxResult;
		}

		approxHappened = transformer.OperationApproximationHappened();

		if (approxHappened || transformer.OperationApproximationHappened())
		{
		    prec *= 4;
		    continue;
		}

		lastBW++;
	    }

	    for (int bw = lastBW; bw <= 32; bw += 2)
	    {
		Result approxResult = runFunction(transformer, bw, prec);
		if (approxResult != UNKNOWN)
		{
		    return approxResult;
		}

		bool approxHappened = transformer.OperationApproximationHappened();
		if (approxHappened || transformer.OperationApproximationHappened())
		{
		    lastBW = bw;
		    break;
		}

		lastBW = bw;
	    }

	    prec *= 4;
	}
    }
    else if (config.approximationMethod == VARIABLES)
    {
	Result approxResult = runFunction(transformer, 1, 0);
	if (approxResult != UNKNOWN)
	{
	    return approxResult;
	}

	approxResult = runFunction(transformer, -1, 0);
	if (approxResult != UNKNOWN)
	{
	    return approxResult;
	}

	for (int bw = 2; bw < 32; bw += 2)
	{
	    approxResult = runFunction(transformer, bw, 0);
	    if (approxResult != UNKNOWN)
	    {
		return approxResult;
	    }
	}
    }

    return UNKNOWN;
}

Result Solver::runWithUnderApproximations(ExprToBDDTransformer &transformer)
{
    return runWithApproximations(transformer, UNDERAPPROXIMATION);
}

Result Solver::runWithOverApproximations(ExprToBDDTransformer &transformer)
{
    return runWithApproximations(transformer, OVERAPPROXIMATION);
}

z3::expr Solver::substituteModel(z3::expr& e, const std::map<std::string, std::vector<bool>>& model) const
{
    auto &context = e.ctx();
    z3::expr_vector consts(context);
    z3::expr_vector vals(context);

    for (auto &varModel : model)
    {
	auto bitwidth = varModel.second.size();
	z3::expr c = context.bv_const(varModel.first.c_str(), bitwidth);

	std::stringstream ss;
	for (auto bit : varModel.second)
	{
	    ss << bit;
	}

	unsigned long long value = 0;
	if (bitwidth <= 8 * sizeof(unsigned long long))
	{
	    value = stoull(ss.str(), 0, 2);
	}

	consts.push_back(c);
	vals.push_back(context.bv_val(value, bitwidth));

	if (bitwidth == 1)
	{
	    consts.push_back(context.bool_const(varModel.first.c_str()));
	    vals.push_back(context.bool_val(varModel.second[0]));
	}
    }

    return e.substitute(consts, vals);
}
