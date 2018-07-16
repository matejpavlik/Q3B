//
// Created by pnavratil on 7/21/17.
//

#ifndef SYLVAN_BVEC_H
#define SYLVAN_BVEC_H

#include <functional>
#include <iostream>
#include <vector>
#include <sylvan_obj.hpp>
#include "../maybeBdd/maybeBdd.h"

namespace sylvan {

class Bvec {
    std::vector<MaybeBDD> m_bitvec;

public:
    Bvec();

    Bvec(size_t bitnum, const Bdd& value);

    Bvec(size_t bitnum, const MaybeBDD& value);

    Bvec(const Bvec& other);

    Bvec& operator=(Bvec other);

    ~Bvec() = default;

    void
    set(size_t i, const Bdd& con);

    void
    set(size_t i, const MaybeBDD& con);

    MaybeBDD&
    operator[](size_t i);

    const MaybeBDD&
    operator[](size_t i) const;

    size_t
    bitnum() const;

    bool
    empty() const;

    static Bvec
    bvec_build(size_t bitnum, bool isTrue);

    static Bvec
    bvec_true(size_t bitnum);

    static Bvec
    bvec_false(size_t bitnum);

    static Bvec
    bvec_con(size_t bitnum, int val);

    static Bvec
    bvec_ncon(size_t bitnum, int val);

    static Bvec
    bvec_var(size_t bitnum, int offset, int step);

    static Bvec
    bvec_varvec(size_t bitnum, int *var);

    Bvec
    bvec_coerce(size_t bitnum) const;

    static Bvec
    arithmetic_neg(const Bvec& src);

    bool
    bvec_isConst() const;

    unsigned int
    bvec_varBits() const;

    int
    bvec_val() const;

    int
    bvec_nval() const;

    void
    bvec_print() const;

    static Bvec
    bvec_copy(const Bvec& other);

    static Bvec
    bvec_map1(const Bvec& src, std::function<MaybeBDD(const MaybeBDD&)> fun);

    static Bvec
    bvec_map2(const Bvec& first, const Bvec& second, std::function<MaybeBDD(const MaybeBDD&, const MaybeBDD&)> fun);

    static Bvec
    bvec_map3(const Bvec& first, const Bvec& second, const Bvec& third, std::function<MaybeBDD(const MaybeBDD&, const MaybeBDD&, const MaybeBDD&)> fun);

    static Bvec
    bvec_add(const Bvec& left, const Bvec& right);

    static Bvec
    bvec_add(const Bvec& left, const Bvec& right, unsigned int);

    static Bvec
    bvec_add_nodeLimit(const Bvec& left, const Bvec& right, unsigned int);

    static Bvec
    bvec_sub(const Bvec& left, const Bvec& right);

    Bvec
    bvec_mulfixed(int con) const;

    static Bvec
    bvec_mul(const Bvec& left, const Bvec& right);

    static Bvec
    bvec_mul(const Bvec& left, const Bvec& right, unsigned int);

    static Bvec
    bvec_mul_nodeLimit(const Bvec& left, const Bvec& right, unsigned int);

    int
    bvec_divfixed(size_t con, Bvec& result, Bvec& rem) const;

    static int
    bvec_div(const Bvec& left, const Bvec& right, Bvec& result, Bvec& rem);

    static int
    bvec_div(const Bvec& left, const Bvec& right, Bvec& result, Bvec& rem, unsigned int);

    static int
    bvec_div_nodeLimit(const Bvec& left, const Bvec& right, Bvec& result, Bvec& rem, unsigned int);

    static Bvec
    bvec_sdiv(const Bvec& left, const Bvec& right);

    static Bvec
    bvec_srem(const Bvec& left, const Bvec& right);

    static Bvec
    bvec_ite(const MaybeBDD& val, const Bvec& left, const Bvec& right);

    static Bvec
    bvec_ite(const MaybeBDD& val, const Bvec& left, const Bvec& right, unsigned int);

    static Bvec
    bvec_ite_nodeLimit(const MaybeBDD& val, const Bvec& left, const Bvec& right, unsigned int);

    Bvec
    bvec_shlfixed(unsigned int pos, const MaybeBDD& con) const;

    static Bvec
    bvec_shl(const Bvec& left, const Bvec& right, const MaybeBDD& con);

    Bvec
    bvec_shrfixed(unsigned int pos, const MaybeBDD& con) const;

    static Bvec
    bvec_shr(const Bvec& left, const Bvec& right, const MaybeBDD& con);

    static MaybeBDD
    bvec_lth(const Bvec& left, const Bvec& right);

    static Bdd
    bvec_lth_overApprox(const Bvec& left, const Bvec& right);

    static Bdd
    bvec_lth_underApprox(const Bvec& left, const Bvec& right);

    static MaybeBDD
    bvec_lte(const Bvec& left, const Bvec& right);

    static Bdd
    bvec_lte_overApprox(const Bvec& left, const Bvec& right);

    static Bdd
    bvec_lte_underApprox(const Bvec& left, const Bvec& right);

    static MaybeBDD
    bvec_gth(const Bvec& left, const Bvec& right);

    static MaybeBDD
    bvec_gte(const Bvec& left, const Bvec& right);

    static MaybeBDD
    bvec_slth(const Bvec& left, const Bvec& right);

    static Bdd
    bvec_slth_overApprox(const Bvec& left, const Bvec& right);

    static Bdd
    bvec_slth_underApprox(const Bvec& left, const Bvec& right);

    static MaybeBDD
    bvec_slte(const Bvec& left, const Bvec& right);

    static Bdd
    bvec_slte_overApprox(const Bvec& left, const Bvec& right);

    static Bdd
    bvec_slte_underApprox(const Bvec& left, const Bvec& right);

    static MaybeBDD
    bvec_sgth(const Bvec& left, const Bvec& right);

    static MaybeBDD
    bvec_sgte(const Bvec& left, const Bvec& right);

    static MaybeBDD
    bvec_equ(const Bvec& left, const Bvec& right);

    static Bdd
    bvec_equ_overApprox(const Bvec& left, const Bvec& right);

    static Bdd
    bvec_equ_underApprox(const Bvec& left, const Bvec& right);

    static MaybeBDD
    bvec_nequ(const Bvec& left, const Bvec& right);

    static Bdd
    bvec_nequ_overApprox(const Bvec& left, const Bvec& right);

    static Bdd
    bvec_nequ_underApprox(const Bvec& left, const Bvec& right);

    Bvec
    operator&(const Bvec& other) const { return bvec_map2(*this, other, bdd_and); }

    Bvec
    operator^(const Bvec& other) const { return bvec_map2(*this, other, bdd_xor); }

    Bvec
    operator|(const Bvec& other) const { return bvec_map2(*this, other, bdd_or); }

    Bvec
    operator!(void) const { return bvec_map1(*this, bdd_not); }

    Bvec
    operator~(void) const { return bvec_map1(*this, bdd_not); }

    Bvec
    operator<<(int con) const { return bvec_shlfixed(con, MaybeBDD(Bdd::bddZero())); }

    Bvec
    operator<<(const Bvec& other) const { return bvec_shl(*this, other, MaybeBDD(Bdd::bddZero())); }

    Bvec
    operator>>(int con) const { return bvec_shrfixed(con, MaybeBDD(Bdd::bddZero())); }

    Bvec
    operator>>(const Bvec& other) const { return bvec_shr(*this, other, MaybeBDD(Bdd::bddZero())); }

    Bvec
    operator+(const Bvec& other) const { return bvec_add(*this, other); }

    Bvec
    operator+=(const Bvec& other) { *this = bvec_add(*this, other); return *this; }


    Bvec
    operator-(const Bvec& other) { return bvec_sub(*this, other); }

    Bvec
    operator-=(const Bvec& other) { *this = bvec_sub(*this, other); return *this; }

    Bvec
    operator*(int con) const { return bvec_mulfixed(con); }

    Bvec
    operator*=(int con) { this->bvec_mulfixed(con); return *this; }

    Bvec
    operator*(const Bvec& other) const { return bvec_mul(*this, other); }

    Bvec
    operator*=(const Bvec& other) { *this = bvec_mul(*this, other); return *this; }

    MaybeBDD
    operator<(const Bvec& other) const { return bvec_lth(*this, other); }

    MaybeBDD
    operator<=(const Bvec& other) const { return bvec_lte(*this, other); }

    MaybeBDD
    operator>(const Bvec& other) const { return bvec_gth(*this, other); }

    MaybeBDD
    operator>=(const Bvec& other) const { return bvec_gte(*this, other); }

    MaybeBDD
    operator==(const Bvec& other) const { return bvec_equ(*this, other); }

    MaybeBDD
    operator!=(const Bvec& other) const { return !(*this == other); }

    unsigned int bddNodes()
    {
	auto count = 0U;

	for (const auto &bdd : m_bitvec)
	{
	    count += bdd.NodeCount();
	}

	return count;
    }

    bool isPrecise() const
    {
	for (const auto &bdd : m_bitvec)
	{
	    if (!bdd.HasValue())
	    {
		return false;
	    }
	}

	return true;
    }

private:

    static void
    bvec_div_rec(Bvec& divisor, Bvec& remainder, Bvec& result, size_t step);

    static MaybeBDD
    bdd_and(const MaybeBDD& first, const MaybeBDD& second);

    static MaybeBDD
    bdd_xor(const MaybeBDD& first, const MaybeBDD& second);

    static MaybeBDD
    bdd_or(const MaybeBDD& first, const MaybeBDD& second);

    static MaybeBDD
    bdd_not(const MaybeBDD& src);

    static MaybeBDD
    get_signs(const MaybeBDD& left, const MaybeBDD& right);

    void
    swap(Bvec& other);

    static Bvec
    reserve(size_t bitnum);

    static void
    reserve(Bvec& bitvector, size_t bitnum);
};

} // sylvan

#endif //BDD_BVEC_H
