#include "maybeBdd.h"

bool MaybeBDD::approximationHappened = false;

MaybeBDD MaybeBDD::And(const MaybeBDD &other) const
{
    if (this->HasValue() && other.HasValue())
    {
	return MaybeBDD(GetBDD() & other.GetBDD());
    }

    if (this->HasValue() && this->GetBDD().isZero())
    {
	return *this;
    }

    if (other.HasValue() && other.GetBDD().isZero())
    {
	return other;
    }

    return MaybeBDD{};
}

MaybeBDD MaybeBDD::Or(const MaybeBDD &other) const
{
    if (this->HasValue() && other.HasValue())
    {
	return MaybeBDD(GetBDD() | other.GetBDD());
    }

    if (this->HasValue() && this->GetBDD().isOne())
    {
	return *this;
    }

    if (other.HasValue() && other.GetBDD().isOne())
    {
	return other;
    }

    return MaybeBDD{};
}

MaybeBDD MaybeBDD::Xor(const MaybeBDD &other) const
{
    if (this->HasValue() && other.HasValue())
    {
	return MaybeBDD(GetBDD() ^ other.GetBDD());
    }

    return MaybeBDD{};
}

MaybeBDD MaybeBDD::Xnor(const MaybeBDD &other) const
{
    if (this->HasValue() && other.HasValue())
    {
	return MaybeBDD(GetBDD().Xnor(other.GetBDD()));
    }

    return MaybeBDD{};
}

MaybeBDD MaybeBDD::Not() const
{
    if (this->HasValue())
    {
	return MaybeBDD(!GetBDD());
    }

    return MaybeBDD{};
}

MaybeBDD MaybeBDD::Ite(const MaybeBDD &thenBdd, const MaybeBDD &elseBdd) const
{
    if (thenBdd.Equals(elseBdd))
    {
	return thenBdd;
    }

    if (this->HasValue() && thenBdd.HasValue() && elseBdd.HasValue())
    {
	return MaybeBDD(this->GetBDD().Ite(thenBdd.GetBDD(), elseBdd.GetBDD()));
    }

    if (this->HasValue() && this->GetBDD().isOne())
    {
	return thenBdd;
    }

    if (this->HasValue() && this->GetBDD().isZero())
    {
	return elseBdd;
    }

    return MaybeBDD{};
}
