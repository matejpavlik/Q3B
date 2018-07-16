#pragma once

#include <sylvan_obj.hpp>

using namespace sylvan;

class BDDInterval
{
public:
    Bdd upper;
    Bdd lower;

    BDDInterval() = default;
    BDDInterval(const BDDInterval&) = default;
    BDDInterval& operator = (const BDDInterval&) = default;

BDDInterval(Bdd bdd) :
    upper(bdd), lower(bdd)
    { }

BDDInterval(Bdd upper, Bdd lower) :
    upper(upper), lower(lower)
    { }

    BDDInterval operator & (const BDDInterval& rhs)
    {
	if (IsPrecise() && rhs.IsPrecise())
	{
	    return BDDInterval{upper & rhs.upper};
	}

	return BDDInterval{upper & rhs.upper, lower & rhs.lower};
    }

    BDDInterval operator * (const BDDInterval& rhs)
    {
	return *this & rhs;
    }

    BDDInterval operator | (const BDDInterval& rhs)
    {
	if (IsPrecise() && rhs.IsPrecise())
	{
	    return BDDInterval{upper | rhs.upper};
	}

	return BDDInterval{upper | rhs.upper, lower | rhs.lower};
    }

    BDDInterval operator + (const BDDInterval& rhs)
    {
	return *this | rhs;
    }

    BDDInterval operator ! ()
    {
	if (IsPrecise())
	{
	    return BDDInterval{!upper};
	}

	return BDDInterval{!upper, !lower};
    }

    BDDInterval Xnor (const BDDInterval& rhs) const
    {
	if (IsPrecise() && rhs.IsPrecise())
	{
	    return BDDInterval{upper.Xnor(rhs.upper)};
	}

	return BDDInterval{upper.Xnor(rhs.upper), lower.Xnor(rhs.lower)};
    }

    BDDInterval Ite (const BDDInterval& t, const BDDInterval& e) const
    {
	if (IsPrecise() && t.IsPrecise() && e.IsPrecise())
	{
	    return BDDInterval{upper.Ite(t.upper, e.upper)};
	}

	return BDDInterval{upper.Ite(t.upper, e.upper), lower.Ite(t.lower, e.lower)};
    }

    BDDInterval UnivAbstract(const BddSet& variables)
    {
	if (IsPrecise())
	{
	    BDDInterval{upper.UnivAbstract(variables)};
	}

	return BDDInterval{upper.UnivAbstract(variables), lower.UnivAbstract(variables)};
    }

    BDDInterval ExistAbstract(const BddSet& variables)
    {
	if (IsPrecise())
	{
	    BDDInterval{upper.ExistAbstract(variables)};
	}

	return BDDInterval{upper.ExistAbstract(variables), lower.ExistAbstract(variables)};
    }

    bool IsPrecise() const
    {
	return upper == lower;
    }
};
