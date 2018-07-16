#pragma once

#include <sylvan_obj.hpp>
#include <optional>
#include <iostream>
#include <vector>

using namespace sylvan;

class MaybeBDD {
private:
    std::optional<Bdd> innerBdd;
    static bool approximationHappened;

public:
MaybeBDD()
    {
	innerBdd = {};
    }

MaybeBDD(Bdd bdd)
    : innerBdd(bdd)
    {

    }

    MaybeBDD& operator=(MaybeBDD other)
    {
	swap(other);
        return *this;
    }

    bool HasValue() const
    {
	return innerBdd.has_value();
    }

    Bdd GetBDD() const
    {
	return innerBdd.value();
    }

    Bdd GetBDD(Bdd ifEmpty) const
    {
	if (!innerBdd.has_value())
	{
	    approximationHappened = true;
	}
	return innerBdd.value_or(ifEmpty);
    }

    unsigned int NodeCount() const
    {
	if (innerBdd.has_value())
	{
	    return innerBdd.value().NodeCount();
	}

	return 0;
    }

    static void ResetApproximationFlag()
    {
	approximationHappened = false;
    }

    static bool ApproximationHappened()
    {
	return approximationHappened;
    }

    MaybeBDD And(const MaybeBDD&) const;
    MaybeBDD Or(const MaybeBDD&) const;
    MaybeBDD Xor(const MaybeBDD&) const;
    MaybeBDD Xnor(const MaybeBDD&) const;
    MaybeBDD Not() const;

    MaybeBDD Ite(const MaybeBDD&, const MaybeBDD&) const;

    bool IsOne() const
    {
	return HasValue() && GetBDD().isOne();
    }

    bool IsZero() const
    {
	return HasValue() && GetBDD().isZero();
    }

    bool IsVar() const
    {
	return !IsZero() && !IsOne();
    }

    MaybeBDD operator&(const MaybeBDD& other) const
    {
	return this->And(other);
    }

    MaybeBDD operator&=(const MaybeBDD& other)
    {
	innerBdd = (*this & other).innerBdd;
	return *this;
    }

    MaybeBDD operator*(const MaybeBDD& other) const
    {
	return this->And(other);
    }

    MaybeBDD operator|(const MaybeBDD& other) const
    {
	return this->Or(other);
    }

    MaybeBDD operator|=(const MaybeBDD& other)
    {
	innerBdd = (*this | other).innerBdd;
	return *this;
    }

    MaybeBDD operator+(const MaybeBDD& other) const
    {
	return this->Or(other);
    }

    MaybeBDD operator!() const
    {
	return this->Not();
    }

    MaybeBDD operator~() const
    {
	return this->Not();
    }

    MaybeBDD operator^(const MaybeBDD& other) const
    {
	return this->Xor(other);
    }

    void swap(MaybeBDD& other)
    {
        using std::swap;
        swap(innerBdd, other.innerBdd);
    }

    operator Bdd() const
    {
	return this->GetBDD();
    }

    bool Equals (const MaybeBDD& other) const
    {
	return (!this->HasValue() && !other.HasValue()) ||
	    (this->HasValue() && other.HasValue() && (this->GetBDD() == other.GetBDD()));
    }
};
