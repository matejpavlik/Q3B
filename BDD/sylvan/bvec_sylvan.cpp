//
// Created by pnavratil on 7/7/17.
//

#include "bvec_sylvan.h"
#include <iostream>

namespace sylvan {

    Bvec::Bvec() : m_bitvec(0) {}

    Bvec::Bvec(size_t bitnum, const MaybeBDD& value)
        : m_bitvec(bitnum, value) {}

    Bvec::Bvec(size_t bitnum, const Bdd& value)
        : m_bitvec(bitnum, MaybeBDD(value)) {}

    Bvec::Bvec(const Bvec& other)
        : m_bitvec(other.m_bitvec) {
    }

    Bvec&
    Bvec::operator=(Bvec other) {
        swap(other);
        return *this;
    }

    void
    Bvec::set(size_t i, const Bdd& con) {
        m_bitvec.at(i) = MaybeBDD(con);
    }

    void
    Bvec::set(size_t i, const MaybeBDD& con) {
        m_bitvec.at(i) = con;
    }

    MaybeBDD&
    Bvec::operator[](size_t i) {
        return m_bitvec.at(i);
    }

    const MaybeBDD&
    Bvec::operator[](size_t i) const {
        return m_bitvec.at(i);
    }

    size_t
    Bvec::bitnum() const {
        return m_bitvec.size();
    }

    bool
    Bvec::empty() const {
        return m_bitvec.empty();
    }

    Bvec
    Bvec::bvec_build(size_t bitnum, bool isTrue) {
        Bvec res(bitnum, isTrue ? Bdd::bddOne() : Bdd::bddZero());
        return res;
    }

    Bvec
    Bvec::bvec_true(size_t bitnum) {
        return bvec_build(bitnum, true);
    }

    Bvec
    Bvec::bvec_false(size_t bitnum) {
        return bvec_build(bitnum, false);
    }

    Bvec
    Bvec::bvec_con(size_t bitnum, int val) {
        Bvec res = reserve(bitnum);
        if (val < 0) {
            throw std::logic_error("use bvec_ncon for negative values");
        }
        for (size_t i = 0U; i < bitnum; ++i) {
            if (val & 1U) {
                res.m_bitvec.push_back(MaybeBDD(Bdd::bddOne()));
            } else {
                res.m_bitvec.push_back(MaybeBDD(Bdd::bddZero()));
            }
            val >>= 1U;
        }
        return res;
    }

    Bvec
    Bvec::bvec_ncon(size_t bitnum, int val) {
        if (val > 0) {
            throw std::logic_error("use bvec_con for positive values");
        }
        val *= -1;
        Bvec res = bvec_con(bitnum, val);
        return arithmetic_neg(res);
    }

    Bvec
    Bvec::bvec_var(size_t bitnum, int offset, int step) {
	LACE_ME;
        Bvec res = bvec_false(bitnum);
        for (size_t i = 0U; i < bitnum; ++i) {
            res[i] = MaybeBDD(Bdd::bddVar(offset + i * step));
        }

        return res;
    }

    Bvec
    Bvec::bvec_varvec(size_t bitnum, int *var) {
        Bvec res = reserve(bitnum);
        for (size_t i = 0U; i < bitnum; ++i) {
            res.m_bitvec.push_back(MaybeBDD(Bdd::bddVar(var[i])));
        }

        return res;
    }

    Bvec
    Bvec::bvec_coerce(size_t bits) const {
        Bvec res = bvec_false(bits);
        size_t minnum = std::min(bits, bitnum());
        for (size_t i = 0U; i < minnum; ++i) {
            res[i] = m_bitvec[i];
        }
        return res;
    }

    Bvec
    Bvec::arithmetic_neg(const Bvec& src) {
        return ~src + bvec_con(src.bitnum(), 1);
    }

    bool
    Bvec::bvec_isConst() const {
        for (size_t i = 0U; i < bitnum(); ++i) {
            if (!(m_bitvec[i].IsOne() || m_bitvec[i].IsZero())) {
                return false;
            }
        }
        return true;
    }

    unsigned int
    Bvec::bvec_varBits() const {
	unsigned int varBits = 0;
        for (size_t i = 0U; i < bitnum(); ++i) {
            if (m_bitvec[i].IsOne() || m_bitvec[i].IsZero()) {
		varBits++;
            }
        }
        return varBits;
    }

    int
    Bvec::bvec_val() const {
        int val = 0;
        for (size_t i = bitnum(); i >= 1U; --i) {
            if (m_bitvec[i - 1U].IsOne())
                val = (val << 1) | 1;
            else if (m_bitvec[i - 1U].IsZero())
                val = val << 1;
            else
                return 0;
        }
        return val;
    }

    int
    Bvec::bvec_nval() const {
        int val = (~(*this) + bvec_con(bitnum(), 1U)).bvec_val();
        return val * -1;
    }

    void
    Bvec::bvec_print() const {
        for (int i = bitnum() - 1; i >= 0; --i) {
	    MaybeBDD currentBdd = m_bitvec[i];
	    if (!currentBdd.HasValue())
	    {
		std::cout << "?";
	    }
	    else
	    {
		if (currentBdd.GetBDD().isZero())
		{
		    std::cout << "0";
		}
		else if (currentBdd.GetBDD().isOne())
		{
		    std::cout << "1";
		}
		else
		{
		    std::cout << "b";
		}
	    }
        }

	std::cout << std::endl;
    }

    Bvec
    Bvec::bvec_copy(const Bvec& other) {
        return Bvec(other);
    }

    Bvec
    Bvec::bvec_map1(const Bvec& src, std::function<MaybeBDD(const MaybeBDD&)> fun) {
        Bvec res = reserve(src.bitnum());
        for (size_t i = 0; i < src.bitnum(); ++i) {
            res.m_bitvec.push_back(fun(src[i]));
        }
        return res;
    }

    Bvec
    Bvec::bvec_map2(const Bvec& first, const Bvec& second, std::function<MaybeBDD(const MaybeBDD&, const MaybeBDD&)> fun) {
	Bvec res;

        if (first.bitnum() != second.bitnum()) {
            return res;
        }

        reserve(res, first.bitnum());
        for (size_t i = 0U; i < first.bitnum(); ++i) {
            res.m_bitvec.push_back(fun(first[i], second[i]));
        }
        return res;
    }

    Bvec
    Bvec::bvec_map3(const Bvec& first, const Bvec& second, const Bvec& third, std::function<MaybeBDD(const MaybeBDD&, const MaybeBDD&, const MaybeBDD&)> fun) {
        Bvec res;

        if (first.bitnum() != second.bitnum() || second.bitnum() != third.bitnum()) {
            return res;
        }

        reserve(res, first.bitnum());
        for (size_t i = 0U; i < first.bitnum(); ++i) {
            res.m_bitvec.push_back(fun(first[i], second[i], third[i]));
        }
        return res;
    }

    Bvec
    Bvec::bvec_add(const Bvec& left, const Bvec& right) {
	return Bvec::bvec_add(left, right, left.bitnum());
    }

    Bvec
    Bvec::bvec_add(const Bvec& left, const Bvec& right, unsigned int precision) {
        Bvec res;
        MaybeBDD comp(Bdd::bddZero());

        if (left.bitnum() == 0 || right.bitnum() == 0 || left.bitnum() != right.bitnum())
        {
            return res;
        }

        reserve(res, left.bitnum());

        for (size_t i = 0U; i < std::min(precision, (unsigned int)left.bitnum()); ++i) {

            /* bitvec[i] = l[i] ^ r[i] ^ c; */
            res.m_bitvec.push_back((left[i] ^ right[i]) ^ comp);

            /* c = (l[i] & r[i]) | (c & (l[i] | r[i])); */
            comp = (left[i] & right[i]) | (comp & (left[i] | right[i]));
        }

	for (size_t i = (size_t)precision; i < left.bitnum(); i++)
	{
	    res.m_bitvec.push_back(MaybeBDD{});
	}

        return res;
    }

    Bvec
    Bvec::bvec_add_nodeLimit(const Bvec& left, const Bvec& right, unsigned int nodeLimit) {
        Bvec res;
        MaybeBDD comp(Bdd::bddZero());

        if (left.bitnum() == 0 || right.bitnum() == 0 || left.bitnum() != right.bitnum())
        {
            return res;
        }

        reserve(res, left.bitnum());

	unsigned int preciseBdds = 0;
        for (size_t i = 0U; i < left.bitnum(); ++i) {

            res.m_bitvec.push_back((left[i] ^ right[i]) ^ comp);

	    preciseBdds++;
	    if (res.bddNodes() > nodeLimit)
	    {
		break;
	    }

            comp = (left[i] & right[i]) | (comp & (left[i] | right[i]));
        }

	for (size_t i = (size_t)preciseBdds; i < left.bitnum(); i++)
	{
	    res.m_bitvec.push_back(MaybeBDD{});
	}

        return res;
    }

    Bvec
    Bvec::bvec_sub(const Bvec& left, const Bvec& right) {
        Bvec res;
        MaybeBDD comp(Bdd::bddZero());

        if (left.bitnum() == 0 || right.bitnum() == 0 || left.bitnum() != right.bitnum())
        {
            return res;
        }

        reserve(res, left.bitnum());

        for (size_t i = 0U; i < left.bitnum(); ++i) {

            /* bitvec[n] = l[n] ^ r[n] ^ comp; */
            res.m_bitvec.push_back((left[i] ^ right[i]) ^ comp);
            /* comp = (l[n] & r[n] & comp) | (!l[n] & (r[n] | comp)); */
            comp = (left[i] & right[i] & comp) | (~left[i] & (right[i] | comp));
        }

        return res;
    }

    Bvec
    Bvec::bvec_mulfixed(int con) const {
        Bvec res;

        if (bitnum() == 0) {
            return res;
        }

        if (con == 0) {
            return bvec_false(bitnum()); /* return false array (base case) */
        }

        Bvec next = bvec_false(bitnum());
        for (size_t i = 1U; i < bitnum(); ++i) {
            next[i] = m_bitvec[i - 1];
        }

        Bvec rest = next.bvec_mulfixed(con >> 1);

        if (con & 0x1) {
            res = bvec_add(*this, rest);
        } else {
            res = rest;
        }

        return res;
    }

    Bvec
    Bvec::bvec_mul(const Bvec& left, const Bvec& right) {
        size_t bitnum = std::max(left.bitnum(), right.bitnum());
	return Bvec::bvec_mul(left, right, bitnum);
    }

    Bvec
    Bvec::bvec_mul(const Bvec& left, const Bvec& right, unsigned int precision)
    {
	size_t bitnum = std::max(left.bitnum(), right.bitnum());
        Bvec res = bvec_false(bitnum);

        if (left.bitnum() == 0 || right.bitnum() == 0) {
            return res;
        }
        Bvec leftshifttmp = Bvec(left);
        Bvec leftshift = leftshifttmp.bvec_coerce(bitnum);

	precision = std::min(precision, (unsigned int)right.bitnum());

        for (size_t i = 0U; i < precision; ++i) {
	    if (!right[i].IsZero())
	    {
		Bvec added = bvec_add(res, leftshift);

		for (size_t m = 0U; m < precision; ++m)
		{
		    res[m] = right[i].Ite(added[m], res[m]);
		}
	    }

            /* Shift 'leftshift' one bit left */
            for (size_t m = bitnum - 1U; m >= 1U; --m) {
                leftshift[m] = leftshift[m - 1];
            }

            //leftshift.m_bitvec.resize(leftshift.bitnum() - 1);
            leftshift[0] = MaybeBDD(Bdd::bddZero());
        }

	for (size_t m = precision; m < bitnum; ++m)
	{
	    res[m] = MaybeBDD{};
	}

        return res;
    }

    Bvec
    Bvec::bvec_mul_nodeLimit(const Bvec& left, const Bvec& right, unsigned int nodeLimit) {
        size_t bitnum = std::max(left.bitnum(), right.bitnum());
        Bvec res = bvec_false(bitnum);

        if (left.bitnum() == 0 || right.bitnum() == 0) {
            return res;
        }
        Bvec leftshifttmp = Bvec(left);
        Bvec leftshift = leftshifttmp.bvec_coerce(bitnum);

	unsigned int preciseBdds = 0;
        for (size_t i = 0U; i < right.bitnum(); ++i) {
	    if (right[i].IsZero())
	    {
		preciseBdds++;
	    }
	    else
	    {
		Bvec added = bvec_add(res, leftshift);

		bool tooManyNodes = false;
		for (size_t m = 0U; m < right.bitnum(); ++m) {

		    res[m] = right[i].Ite(added[m], res[m]);

		    if (res[m].NodeCount() > nodeLimit)
		    {
			tooManyNodes = true;

			if (m >= preciseBdds)
			{
			    preciseBdds++;
			}

			break;
		    }
		}

		if (tooManyNodes)
		{
		    break;
		}
		else
		{
		    preciseBdds++;
		}
	    }

            /* Shift 'leftshift' one bit left */
            for (size_t m = bitnum - 1U; m >= 1U; --m) {
                leftshift[m] = leftshift[m - 1];
            }

            leftshift[0] = MaybeBDD(Bdd::bddZero());
        }

	for (size_t m = preciseBdds; m < bitnum; ++m)
	{
	    res[m] = MaybeBDD{};
	}

        return res;
    }

    void
    Bvec::bvec_div_rec(Bvec& divisor, Bvec& remainder, Bvec& result, size_t step) {
        MaybeBDD isSmaller = bvec_lte(divisor, remainder);
        Bvec newResult = result.bvec_shlfixed(1, isSmaller);
        Bvec zero = bvec_build(divisor.bitnum(), false);
        Bvec sub = reserve(divisor.bitnum());

        for (size_t i = 0U; i < divisor.bitnum(); ++i) {
            sub.m_bitvec.push_back(isSmaller.Ite(divisor[i], zero[i]));
        }

        Bvec tmp = remainder - sub;
        Bvec newRemainder = tmp.bvec_shlfixed(1, result[divisor.bitnum() - 1]);

        if (step > 1) {
            bvec_div_rec(divisor, newRemainder, newResult, step - 1);
        }

        result = newResult;
        remainder = newRemainder;
    }

    int
    Bvec::bvec_divfixed(size_t con, Bvec& result, Bvec& rem) const {
        if (con > 0) {
            Bvec divisor = bvec_con(bitnum(), con);
            Bvec tmp = bvec_false(bitnum());
            Bvec tmpremainder = tmp.bvec_shlfixed(1, m_bitvec[bitnum() - 1]);
            Bvec res = bvec_shlfixed(1, MaybeBDD(Bdd::bddZero()));

            bvec_div_rec(divisor, tmpremainder, result, divisor.bitnum());
            Bvec remainder = tmpremainder.bvec_shrfixed(1, MaybeBDD(Bdd::bddZero()));

            result = res;
            rem = remainder;
            return 0;
        }
        return 1;
    }

    int
    Bvec::bvec_div(const Bvec& left, const Bvec& right, Bvec& result, Bvec& remainder) {
        size_t bitnum = left.bitnum() + right.bitnum();
	return Bvec::bvec_div(left, right, result, remainder, bitnum);
    }

    int
    Bvec::bvec_div(const Bvec& left, const Bvec& right, Bvec& result, Bvec& remainder, unsigned int precision) {
        size_t bitnum = left.bitnum() + right.bitnum();
        if (left.bitnum() == 0 || right.bitnum() == 0 || left.bitnum() != right.bitnum()) {
            return 1;
        }

        Bvec rem = left.bvec_coerce(bitnum);
        Bvec divtmp = right.bvec_coerce(bitnum);
        Bvec div = divtmp.bvec_shlfixed(left.bitnum(), MaybeBDD(Bdd::bddZero()));

        Bvec res = bvec_false(right.bitnum());

        for (size_t i = 0U; i < right.bitnum() + 1 && i < precision; ++i)
	{
            MaybeBDD divLteRem = bvec_lte(div, rem);
            Bvec remSubDiv = bvec_sub(rem, div);

            for (size_t j = 0U; j < bitnum; ++j) {
                rem[j] = divLteRem.Ite(remSubDiv[j], rem[j]);
            }

            if (i > 0) {
                res[right.bitnum() - i] = divLteRem;
            }

            /* Shift 'div' one bit right */
            for (size_t j = 0U; j < bitnum - 1; ++j) {
                div[j] = div[j + 1];
            }
            div[bitnum - 1] = MaybeBDD(Bdd::bddZero());
        }

	//forget lower bits, as then can be imprecise
	for (unsigned int i = precision; i < right.bitnum(); i++)
	{
	    res[right.bitnum() - i - 1] = MaybeBDD{};
	    rem[right.bitnum() - i - 1] = MaybeBDD{};
	}

        result = res;
        remainder = rem.bvec_coerce(right.bitnum());
        return 0;
    }

    int
    Bvec::bvec_div_nodeLimit(const Bvec& left, const Bvec& right, Bvec& result, Bvec& remainder, unsigned int nodeLimit) {
        size_t bitnum = left.bitnum() + right.bitnum();
        if (left.bitnum() == 0 || right.bitnum() == 0 || left.bitnum() != right.bitnum()) {
            return 1;
        }

        Bvec rem = left.bvec_coerce(bitnum);
        Bvec divtmp = right.bvec_coerce(bitnum);
        Bvec div = divtmp.bvec_shlfixed(left.bitnum(), MaybeBDD(Bdd::bddZero()));

        Bvec res = bvec_false(right.bitnum());

	unsigned int preciseBdds = 0;
        for (size_t i = 0U; i < right.bitnum() + 1; ++i)
	{
            MaybeBDD divLteRem = bvec_lte(div, rem);
            Bvec remSubDiv = bvec_sub(rem, div);

            for (size_t j = 0U; j < bitnum; ++j) {
                rem[j] = divLteRem.Ite(remSubDiv[j], rem[j]);
            }

            if (i > 0) {
                res[right.bitnum() - i] = divLteRem;
            }

	    preciseBdds++;
	    if (res.bddNodes() > nodeLimit || rem.bddNodes() > nodeLimit)
	    {
		break;
	    }

	    /* Shift 'div' one bit right */
            for (size_t j = 0U; j < bitnum - 1; ++j) {
                div[j] = div[j + 1];
            }

            div[bitnum - 1] = MaybeBDD(Bdd::bddZero());
        }

	//the first bit of the result was not stored
	if (preciseBdds > 0)
	{
	    preciseBdds--;
	}

	//forget lower bits, as then can be imprecise
	for (unsigned int i = preciseBdds; i < right.bitnum(); i++)
	{
	    res[right.bitnum() - i - 1] = MaybeBDD{};
	}

	if (preciseBdds != right.bitnum())
	{
	    for (unsigned int i = 0; i < right.bitnum(); i++)
	    {
		rem[i] = MaybeBDD{};
	    }
	}

        result = res;
        remainder = rem.bvec_coerce(right.bitnum());
        return 0;
    }

    Bvec
    Bvec::bvec_sdiv(const Bvec& left, const Bvec& right) {
        if (left.bitnum() == 0 || right.bitnum() == 0 || left.bitnum() != right.bitnum()) {
	    abort();
        }
        size_t size = left.bitnum() - 1;

	const MaybeBDD& lhead = left[size];
        const MaybeBDD& rhead = right[size];

	Bvec nnDiv = bvec_false(left.bitnum());
	Bvec nnRem = nnDiv;
	Bvec pnDiv = nnDiv;
	Bvec pnRem = nnDiv;
	Bvec npDiv = nnDiv;
	Bvec npRem = nnDiv;
	Bvec ppDiv = nnDiv;
	Bvec ppRem = nnDiv;

	bvec_div(left, right, nnDiv, nnRem);
	bvec_div(arithmetic_neg(left), right, pnDiv, pnRem);
	bvec_div(left, arithmetic_neg(right), npDiv, npRem);
	bvec_div(arithmetic_neg(left), arithmetic_neg(right), ppDiv, ppRem);

	return bvec_ite((!lhead) & (!rhead),
			nnDiv,
		        bvec_ite(lhead & !rhead,
				 arithmetic_neg(pnDiv),
				 bvec_ite((!lhead) & rhead,
					  arithmetic_neg(npDiv),
					  ppDiv)));
    }

    Bvec
    Bvec::bvec_srem(const Bvec& left, const Bvec& right) {
        if (left.bitnum() == 0 || right.bitnum() == 0 || left.bitnum() != right.bitnum()) {
	    abort();
        }
        size_t size = left.bitnum() - 1;

	const MaybeBDD& lhead = left[size];
        const MaybeBDD& rhead = right[size];

	Bvec nnDiv = bvec_false(left.bitnum());
	Bvec nnRem = nnDiv;
	Bvec pnDiv = nnDiv;
	Bvec pnRem = nnDiv;
	Bvec npDiv = nnDiv;
	Bvec npRem = nnDiv;
	Bvec ppDiv = nnDiv;
	Bvec ppRem = nnDiv;

	bvec_div(left, right, nnDiv, nnRem);
	bvec_div(arithmetic_neg(left), right, pnDiv, pnRem);
	bvec_div(left, arithmetic_neg(right), npDiv, npRem);
	bvec_div(arithmetic_neg(left), arithmetic_neg(right), ppDiv, ppRem);

	return bvec_ite((!lhead) & !rhead,
			nnRem,
		        bvec_ite(lhead & !rhead,
				 arithmetic_neg(pnRem),
				 bvec_ite((!lhead) & rhead,
					  npRem,
					  arithmetic_neg(ppRem))));
    }

    Bvec
    Bvec::bvec_ite(const MaybeBDD& val, const Bvec& left, const Bvec& right) {
	return Bvec::bvec_ite(val, left, right, left.bitnum());
    }

    Bvec
    Bvec::bvec_ite(const MaybeBDD& val, const Bvec& left, const Bvec& right, unsigned int precision) {
        Bvec res;

        if (left.bitnum() != right.bitnum()) {
            return res;
        }
        reserve(res, left.bitnum());
        for (size_t i = 0U; i < left.bitnum() && i < precision; ++i) {
            res.m_bitvec.push_back(val.Ite(left[i], right[i]));
        }

	for (size_t i = precision; i < left.bitnum(); ++i) {
	    res.m_bitvec.push_back(MaybeBDD{});
        }
        return res;
    }

    Bvec
    Bvec::bvec_ite_nodeLimit(const MaybeBDD& val, const Bvec& left, const Bvec& right, unsigned int nodeLimit) {
        Bvec res;

        if (left.bitnum() != right.bitnum()) {
            return res;
        }
        reserve(res, left.bitnum());

	auto preciseBdds = 0U;
	if (nodeLimit != 0)
	{
	    for (size_t i = 0U; i < left.bitnum(); ++i) {
		res.m_bitvec.push_back(val.Ite(left[i], right[i]));

		preciseBdds++;

		if (res.bddNodes() > nodeLimit)
		{
		    break;
		}
	    }
	}

	for (size_t i = preciseBdds; i < left.bitnum(); ++i) {
	    res.m_bitvec.push_back(MaybeBDD{});
        }
        return res;
    }

    Bvec
    Bvec::bvec_shlfixed(unsigned int pos, const MaybeBDD& con) const {

        size_t min = (bitnum() < pos) ? bitnum() : pos;

        if (pos < 0U || bitnum() == 0) {
            return Bvec{};
        }

        Bvec res(bitnum(), con);
        for (size_t i = min; i < bitnum(); i++) {
            res[i] = m_bitvec[i - pos];
        }
        return res;
    }

    Bvec
    Bvec::bvec_shl(const Bvec& left, const Bvec& right, const MaybeBDD& con) {
        Bvec res;
        if (left.bitnum() == 0 || right.bitnum() == 0) {

            return res;
        }

        res = bvec_false(left.bitnum());

        for (size_t i = 0U; i <= left.bitnum(); ++i) {
            Bvec val = bvec_con(right.bitnum(), i);
            MaybeBDD rEquN = bvec_equ(right, val);

            for (size_t j = 0U; j < left.bitnum(); ++j) {
                /* Set the m'th new location to be the (m+n)'th old location */
                if (j >= i) {
                    res[j] |= rEquN & left[j - i];
                } else {
		    res[j] |= rEquN & con;
                }
	    }
        }

        /* At last make sure 'c' is shifted in for r-values > l-bitnum */
        Bvec val = bvec_con(right.bitnum(), left.bitnum());
        MaybeBDD rEquN = bvec_gth(right, val);

        for (size_t i = 0U; i < left.bitnum(); i++) {
            res[i] |= (rEquN & con);
        }

        return res;
    }

    Bvec
    Bvec::bvec_shrfixed(unsigned int pos, const MaybeBDD& con) const {
        if (pos < 0 || bitnum() == 0) {
            return Bvec{};
        }
        unsigned int maxnum = std::max(static_cast<unsigned int>(bitnum()) - pos, 0U);
        Bvec res(bitnum(), con);

        for (size_t i = 0U; i < maxnum; ++i) {
            res[i] = m_bitvec[i + pos];
        }
        return res;
    }

    Bvec
    Bvec::bvec_shr(const Bvec& left, const Bvec& right, const MaybeBDD& con) {
        Bvec res;
        if (left.bitnum() == 0 || right.bitnum() == 0) {
            return res;
        }

        res = bvec_false(left.bitnum());
        MaybeBDD tmp1, rEquN;

        for (size_t i = 0U; i <= left.bitnum(); ++i) {
            Bvec val = bvec_con(right.bitnum(), i);
            rEquN = right == val;

            for (size_t j = 0U; j < left.bitnum(); ++j) {
                /* Set the m'th new location to be the (m+n)'th old location */
                if (j + i < left.bitnum())
                    tmp1 = rEquN & left[j + i];
                else
                    tmp1 = rEquN & con;
                res[j] = res[j] | tmp1;
            }
        }

        /* At last make sure 'c' is shifted in for r-values > l-bitnum */
        Bvec val = bvec_con(right.bitnum(), left.bitnum());
        rEquN = bvec_gth(right, val);
        tmp1 = rEquN & con;

        for (size_t i = 0U; i < left.bitnum(); ++i) {
            res[i] = res[i] | tmp1;
        }
        return res;
    }

    MaybeBDD
    Bvec::bvec_lth(const Bvec& left, const Bvec& right) {
        MaybeBDD p(Bdd::bddZero());

        if (left.bitnum() == 0 || right.bitnum() == 0 || left.bitnum() != right.bitnum()) {
            return p;
        }

        for (size_t i = 0U; i < left.bitnum(); ++i) {
            /* p = (!l[n] & r[n]) |
             *     bdd_apply(l[n], r[n], bddop_biimp) & p; */
            p = ((!left[i]) & right[i]) | (left[i].Xnor(right[i]) & p);
        }

        return p;
    }

    Bdd
    Bvec::bvec_lth_overApprox(const Bvec& left, const Bvec& right) {
        Bdd p = Bdd::bddZero();

        if (left.bitnum() == 0 || right.bitnum() == 0 || left.bitnum() != right.bitnum()) {
            return p;
        }

        for (size_t i = 0U; i < left.bitnum(); ++i) {
            /* p = (!l[n] & r[n]) |
             *     bdd_apply(l[n], r[n], bddop_biimp) & p; */
            p = ((!left[i]) & right[i]).GetBDD(Bdd::bddOne()) |
                (left[i].Xnor(right[i]).GetBDD(Bdd::bddOne()) & p);
        }

        return p;
    }

    Bdd
    Bvec::bvec_lth_underApprox(const Bvec& left, const Bvec& right) {
        Bdd p = Bdd::bddZero();

        if (left.bitnum() == 0 || right.bitnum() == 0 || left.bitnum() != right.bitnum()) {
            return p;
        }

        for (size_t i = 0U; i < left.bitnum(); ++i) {
            /* p = (!l[n] & r[n]) |
             *     bdd_apply(l[n], r[n], bddop_biimp) & p; */
            p = ((!left[i]) & right[i]).GetBDD(Bdd::bddZero()) |
                (left[i].Xnor(right[i]).GetBDD(Bdd::bddZero()) & p);
        }

        return p;
    }

    MaybeBDD
    Bvec::bvec_lte(const Bvec& left, const Bvec& right) {
        MaybeBDD p(Bdd::bddOne());

        if (left.bitnum() == 0 || right.bitnum() == 0 || left.bitnum() != right.bitnum()) {
            return p;
        }

        for (size_t i = 0U; i < left.bitnum(); ++i) {
            /* p = (!l[n] & r[n]) |
             *     bdd_apply(l[n], r[n], bddop_biimp) & p; */
            p = ((!left[i]) & right[i]) | (left[i].Xnor(right[i]) & p);
        }

        return p;
    }

    Bdd
    Bvec::bvec_lte_overApprox(const Bvec& left, const Bvec& right) {
        Bdd p = Bdd::bddOne();

        if (left.bitnum() == 0 || right.bitnum() == 0 || left.bitnum() != right.bitnum()) {
            return p;
        }

        for (size_t i = 0U; i < left.bitnum(); ++i) {
            /* p = (!l[n] & r[n]) |
             *     bdd_apply(l[n], r[n], bddop_biimp) & p; */
            p = ((!left[i]) & right[i]).GetBDD(Bdd::bddOne()) |
                (left[i].Xnor(right[i]).GetBDD(Bdd::bddOne()) & p);
        }

        return p;
    }

    Bdd
    Bvec::bvec_lte_underApprox(const Bvec& left, const Bvec& right) {
        Bdd p = Bdd::bddOne();

        if (left.bitnum() == 0 || right.bitnum() == 0 || left.bitnum() != right.bitnum()) {
            return p;
        }

        for (size_t i = 0U; i < left.bitnum(); ++i) {
            /* p = (!l[n] & r[n]) |
             *     bdd_apply(l[n], r[n], bddop_biimp) & p; */
            p = ((!left[i]) & right[i]).GetBDD(Bdd::bddZero()) |
                (left[i].Xnor(right[i]).GetBDD(Bdd::bddZero()) & p);
        }

        return p;
    }

    MaybeBDD
    Bvec::bvec_gth(const Bvec& left, const Bvec& right) {
        return bvec_lth(right, left);
    }

    MaybeBDD
    Bvec::bvec_gte(const Bvec& left, const Bvec& right) {
        return !bvec_lte(right, left);
    }

    MaybeBDD
    Bvec::bvec_slth(const Bvec& left, const Bvec& right) {

        if (left.bitnum() == 0 || right.bitnum() == 0) {
            return MaybeBDD(Bdd::bddZero());
        }

        size_t size = left.bitnum() - 1;

        return get_signs(left[size], right[size]) |
	    (left[size].Xnor(right[size]) &
	     bvec_lth(left.bvec_coerce(size), right.bvec_coerce(size)));
    }

    Bdd
    Bvec::bvec_slth_overApprox(const Bvec& left, const Bvec& right) {

        if (left.bitnum() == 0 || right.bitnum() == 0) {
            return Bdd::bddZero();
        }

        size_t size = left.bitnum() - 1;

        return get_signs(left[size], right[size]).GetBDD(Bdd::bddOne()) |
	    (left[size].Xnor(right[size]).GetBDD(Bdd::bddOne()) &
	     bvec_lth_overApprox(left.bvec_coerce(size), right.bvec_coerce(size)));
    }


    Bdd
    Bvec::bvec_slth_underApprox(const Bvec& left, const Bvec& right) {

        if (left.bitnum() == 0 || right.bitnum() == 0) {
            return Bdd::bddZero();
        }

        size_t size = left.bitnum() - 1;

        return get_signs(left[size], right[size]).GetBDD(Bdd::bddZero()) |
	    (left[size].Xnor(right[size]).GetBDD(Bdd::bddZero()) &
	     bvec_lth_underApprox(left.bvec_coerce(size), right.bvec_coerce(size)));
    }

    MaybeBDD
    Bvec::bvec_slte(const Bvec& left, const Bvec& right) {
        if (left.bitnum() == 0 || right.bitnum() == 0) {
            return MaybeBDD(Bdd::bddZero());
        }

        size_t size = left.bitnum() - 1;

        return (get_signs(left[size], right[size]) |
		(left[size].Xnor(right[size]) &
		 bvec_lte(left.bvec_coerce(size), right.bvec_coerce(size))));
    }

    Bdd
    Bvec::bvec_slte_overApprox(const Bvec& left, const Bvec& right) {
        if (left.bitnum() == 0 || right.bitnum() == 0) {
            return Bdd::bddZero();
        }

        size_t size = left.bitnum() - 1;

        return get_signs(left[size], right[size]).GetBDD(Bdd::bddOne()) |
	    (left[size].Xnor(right[size]).GetBDD(Bdd::bddOne()) &
	     bvec_lte_overApprox(left.bvec_coerce(size), right.bvec_coerce(size)));
    }

    Bdd
    Bvec::bvec_slte_underApprox(const Bvec& left, const Bvec& right) {
        if (left.bitnum() == 0 || right.bitnum() == 0) {
            return Bdd::bddZero();
        }

        size_t size = left.bitnum() - 1;

        return get_signs(left[size], right[size]).GetBDD(Bdd::bddZero()) |
	    (left[size].Xnor(right[size]).GetBDD(Bdd::bddZero()) &
	     bvec_lte_underApprox(left.bvec_coerce(size), right.bvec_coerce(size)));
    }

    MaybeBDD
    Bvec::get_signs(const MaybeBDD& left, const MaybeBDD& right) {
        MaybeBDD differentSigns =
	    left.Xnor(MaybeBDD(Bdd::bddOne())) &
	    right.Xnor(MaybeBDD(Bdd::bddZero()));
        return differentSigns;
    }

    MaybeBDD
    Bvec::bvec_sgth(const Bvec& left, const Bvec& right) {
        return !bvec_slte(left, right);
    }

    MaybeBDD
    Bvec::bvec_sgte(const Bvec& left, const Bvec& right) {
        return !bvec_slth(left, right);
    }

    MaybeBDD
    Bvec::bvec_equ(const Bvec& left, const Bvec& right) {
        MaybeBDD p(Bdd::bddOne());

        if (left.bitnum() == 0 || right.bitnum() == 0 || left.bitnum() != right.bitnum()) {
            return MaybeBDD(Bdd::bddZero());
        }

        for (size_t i = 0U; i < left.bitnum(); ++i) {
            p = p & left[i].Xnor(right[i]);
	    if (p.IsZero())
	    {
		return p;
	    }
        }
        return p;
    }

    Bdd
    Bvec::bvec_equ_overApprox(const Bvec& left, const Bvec& right) {
        Bdd p = Bdd::bddOne();

        if (left.bitnum() == 0 || right.bitnum() == 0 || left.bitnum() != right.bitnum()) {
            return Bdd::bddZero();
        }

        for (size_t i = 0U; i < left.bitnum(); ++i) {
            p = p & left[i].Xnor(right[i]).GetBDD(Bdd::bddOne());
	    if (p.isZero())
	    {
		return p;
	    }
        }
        return p;
    }

    Bdd
    Bvec::bvec_equ_underApprox(const Bvec& left, const Bvec& right) {
        Bdd p = Bdd::bddOne();

        if (left.bitnum() == 0 || right.bitnum() == 0 || left.bitnum() != right.bitnum()) {
            return Bdd::bddZero();
        }

        for (size_t i = 0U; i < left.bitnum(); ++i) {
            p = p & left[i].Xnor(right[i]).GetBDD(Bdd::bddZero());
	    if (p.isZero())
	    {
		return p;
	    }
        }
        return p;
    }

    MaybeBDD
    Bvec::bvec_nequ(const Bvec& left, const Bvec& right) {
        MaybeBDD p(Bdd::bddZero());

        if (left.bitnum() == 0 || right.bitnum() == 0 || left.bitnum() != right.bitnum()) {
            return MaybeBDD(Bdd::bddZero());
        }

        for (size_t i = 0U; i < left.bitnum(); ++i) {
            p = p | left[i].Xor(right[i]);
	    if (p.IsOne())
	    {
		return p;
	    }
        }
        return p;
    }

    Bdd
    Bvec::bvec_nequ_overApprox(const Bvec& left, const Bvec& right) {
        Bdd p = Bdd::bddZero();

        if (left.bitnum() == 0 || right.bitnum() == 0 || left.bitnum() != right.bitnum()) {
            return Bdd::bddZero();
        }

        for (size_t i = 0U; i < left.bitnum(); ++i) {
            p = p | left[i].Xor(right[i]).GetBDD(Bdd::bddOne());
	    if (p.isOne())
	    {
		return p;
	    }
        }
        return p;
    }

    Bdd
    Bvec::bvec_nequ_underApprox(const Bvec& left, const Bvec& right) {
        Bdd p = Bdd::bddZero();

        if (left.bitnum() == 0 || right.bitnum() == 0 || left.bitnum() != right.bitnum()) {
            return Bdd::bddZero();
        }

        for (size_t i = 0U; i < left.bitnum(); ++i) {
            p = p | (left[i].Xor(right[i])).GetBDD(Bdd::bddZero());
	    if (p.isOne())
	    {
		return p;
	    }
        }
        return p;
    }

    MaybeBDD
    Bvec::bdd_and(const MaybeBDD& first, const MaybeBDD& second) {
        return first & second;
    }

    MaybeBDD
    Bvec::bdd_xor(const MaybeBDD& first, const MaybeBDD& second) {
        return first ^ second;
    }

    MaybeBDD
    Bvec::bdd_or(const MaybeBDD& first, const MaybeBDD& second) {
        return first | second;
    }

    MaybeBDD
    Bvec::bdd_not(const MaybeBDD& src) {
        return !src;
    }

    void
    Bvec::swap(Bvec& other) {
        using std::swap;
        swap(m_bitvec, other.m_bitvec);
    }

    Bvec
    Bvec::reserve(size_t bitnum) {
        Bvec res;
        res.m_bitvec.reserve(bitnum);
        return res;
    }

    void
    Bvec::reserve(Bvec& bitvector, size_t bitnum) {
        bitvector.m_bitvec.reserve(bitnum);
    }

} //sylvan
