#ifndef BDD_BVEC_H
#define BDD_BVEC_H

#include <chrono>
#include <assert.h>
#include <functional>
#include <fstream>
#include <vector>
#include <cuddObj.hh>
#include <iostream>

using namespace std::chrono;

namespace cudd {

class Bvec {
    Cudd* m_manager;

public:
    std::vector<BDD> m_bitvec;

    Bvec() = delete;

    Bvec(Cudd& manager);

    Bvec(Cudd& manager, size_t bitnum, const BDD& value);

    Bvec(const Bvec& other);

    Bvec& operator=(Bvec other);

    ~Bvec() = default;

    void
    set(size_t i, const BDD& con);

    BDD&
    operator[](size_t i);

    const BDD&
    operator[](size_t i) const;

    size_t
    bitnum() const;

    Cudd&
    manager() const;

    bool
    empty() const;

    static Bvec
    bvec_build(Cudd& manager, size_t bitnum, bool isTrue);

    static Bvec
    bvec_true(Cudd& manager, size_t bitnum);

    static Bvec
    bvec_false(Cudd& manager, size_t bitnum);

    static Bvec
    bvec_con(Cudd& manager, size_t bitnum, int val);

    static Bvec
    bvec_var(Cudd& manager, size_t bitnum, int offset, int step);

    Bvec
    bvec_coerce(size_t bitnum) const;

    static Bvec
    arithmetic_neg(const Bvec& src);

    bool
    bvec_isPreciseConst() const;

    unsigned int
    bvec_varBits() const;

    int
    bvec_val() const;

    static Bvec
    bvec_copy(const Bvec& other);

    static Bvec
    bvec_map1(const Bvec& src, const std::function<BDD(const BDD&)>& fun);

    static Bvec
    bvec_map2(const Bvec& first, const Bvec& second, const std::function<BDD(const BDD&, const BDD&)>& fun);

    static Bvec
    bvec_add(const Bvec& left, const Bvec& right, bool precise);

    static Bvec
    bvec_add_nodeLimit(const Bvec& left, const Bvec& right, bool precise, unsigned int);

    static Bvec
    bvec_add_reduced(const Bvec& left, const Bvec& right, traverse_heuristic heu, unsigned int nodeLimit);

    static Bvec
    bvec_sub(const Bvec& left, const Bvec& right, bool precise);

    static Bvec
    bvec_sub_reduced(const Bvec& left, const Bvec& right, traverse_heuristic heu, unsigned int nodeLimit);

    Bvec
    bvec_mulfixed(int con, bool precise) const;

    static Bvec
    bvec_mul(const Bvec& left, const Bvec& right, bool precise);

    static Bvec
    bvec_mul_nodeLimit(const Bvec& left, const Bvec& right, bool precise, unsigned int);

    static Bvec
    bvec_mul_reduced(const Bvec& left, const Bvec& right, traverse_heuristic heu, unsigned int nodeLimit);

    int
    bvec_divfixed(size_t con, Bvec& result, Bvec& rem, bool precise) const;

    static int
    bvec_div(const Bvec& left, const Bvec& right, Bvec& result, Bvec& rem, bool precise);

    static int
    bvec_div_nodeLimit(const Bvec& left, const Bvec& right, Bvec& result, Bvec& rem, bool precise, unsigned int);

    static int
    bvec_div_reduced(const Bvec& left, const Bvec& right, Bvec& result, Bvec& rem, traverse_heuristic heu, unsigned int nodeLimit);

    static Bvec
    bvec_ite(const BDD& val, const Bvec& left, const Bvec& right, bool precise);

    static Bvec
    bvec_ite_nodeLimit(const BDD& val, const Bvec& left, const Bvec& right, bool precise, unsigned int);

    static Bvec
    bvec_ite_reduced(const BDD& val, const Bvec& left, const Bvec& right, traverse_heuristic heu, unsigned int nodeLimit);

    Bvec
    bvec_shlfixed(unsigned int pos, const BDD& con) const;

    static Bvec
    bvec_shl(const Bvec& left, const Bvec& right, const BDD& con, bool precise);

    static Bvec
    bvec_shl_reduced(const Bvec& left, const Bvec& right, const BDD& con, traverse_heuristic heu, unsigned int nodeLimit);

    Bvec
    bvec_shrfixed(unsigned int pos, const BDD& con) const;

    static Bvec
    bvec_shr(const Bvec& left, const Bvec& right, const BDD& con, bool precise);

    static Bvec
    bvec_shr_reduced(const Bvec& left, const Bvec& right, const BDD& con, traverse_heuristic heu, unsigned int nodeLimit);

    static BDD
    bvec_lth(const Bvec& left, const Bvec& right, bool precise) {
        Cudd& manager = check_same_cudd(*left.m_manager, *right.m_manager);
        BDD p =  manager.bddZero();

        if (left.bitnum() == 0 || right.bitnum() == 0 || left.bitnum() != right.bitnum()) {
            return p;
        }
        
        for (size_t i = 0U; i < left.bitnum(); ++i) {
            if (precise) {
                p = ((~left[i]).AndP(right[i])).OrP(left[i].XnorP(right[i]).AndP(p));
            } else {
                p = (~left[i] & right[i]) | (left[i].Xnor(right[i]) & p);
                
                /* alternative 1 
                if (!right[i].IsUnknown()) {
                   p = right[i].Ite(~left[i] | p, ~left[i] & p);
                } else if (!left[i].IsUnknown()) {
                   p = left[i].Ite(right[i] & p, right[i] | p);
                } else {
                   p = manager.bddUnknown();
                }*/
                
                /* alternative 2 
                if (right[i].IsOne()) {
                   p |= ~left[i];
                } else if (right[i].IsZero()) {
                   p &= ~left[i];
                } else if (left[i].IsOne()) {
                   p &= right[i];
                } else if (left[i].IsZero()) {
                   p |= right[i];
                } else {
                   p = (~left[i] & right[i]) | (left[i].Xnor(right[i]) & p);
                }
                */
            }
        }

        return p;
    }

    static BDD
    bvec_lth_reduced(const Bvec& left, const Bvec& right, traverse_heuristic heu, unsigned int nodeLimit) {
        Cudd& manager = check_same_cudd(*left.m_manager, *right.m_manager);
        BDD p =  manager.bddZero();

        if (left.bitnum() == 0 || right.bitnum() == 0 || left.bitnum() != right.bitnum()) {
            return p;
        }
        
        for (size_t i = 0U; i < left.bitnum(); ++i) {
            p = ((~left[i]).AndLim(right[i], heu, nodeLimit))
                .OrLim(left[i].XnorLim(right[i], heu, nodeLimit).AndLim(p, heu, nodeLimit), heu, nodeLimit);
        }

        return p;
    }

    static BDD
    bvec_lte(const Bvec& left, const Bvec& right, bool precise) {
        Cudd& manager = check_same_cudd(*left.m_manager, *right.m_manager);
        BDD p = manager.bddOne();

        if (left.bitnum() == 0 || right.bitnum() == 0 || left.bitnum() != right.bitnum()) {
            return p;
        }
        
        for (size_t i = 0U; i < left.bitnum(); ++i) {
            if (precise) {
                p = ((~left[i]).AndP(right[i])).OrP(left[i].XnorP(right[i]).AndP(p));
            } else {
                p = (~left[i] & right[i]) | (left[i].Xnor(right[i]) & p);
                
                /* alternative 1 
                if (!right[i].IsUnknown()) {
                   p = right[i].Ite(~left[i] | p, ~left[i] & p);
                } else if (!left[i].IsUnknown()) {
                   p = left[i].Ite(right[i] & p, right[i] | p);
                } else {
                   p = manager.bddUnknown();
                } */
                
                
                /* alternative 2 
                if (right[i].IsOne()) {
                   p |= ~left[i];
                } else if (right[i].IsZero()) {
                   p &= ~left[i];
                } else if (left[i].IsOne()) {
                   p &= right[i];
                } else if (left[i].IsZero()) {
                   p |= right[i];
                } else {
                   p = (~left[i] & right[i]) | (left[i].Xnor(right[i]) & p);
                }
                */
            }
        }

        return p;
    }

    static BDD
    bvec_lte_reduced(const Bvec& left, const Bvec& right, traverse_heuristic heu, unsigned int nodeLimit) {
        Cudd& manager = check_same_cudd(*left.m_manager, *right.m_manager);
        BDD p = manager.bddOne();

        if (left.bitnum() == 0 || right.bitnum() == 0 || left.bitnum() != right.bitnum()) {
            return p;
        }
        
        for (size_t i = 0U; i < left.bitnum(); ++i) {
            p = ((~left[i]).AndLim(right[i], heu, nodeLimit))
                .OrLim(left[i].XnorLim(right[i], heu, nodeLimit).AndLim(p, heu, nodeLimit), heu, nodeLimit);
        }

        return p;
    }

    static BDD
    bvec_gth(const Bvec& left, const Bvec& right, bool precise);

    static BDD
    bvec_gth_reduced(const Bvec& left, const Bvec& right, traverse_heuristic heu, unsigned int nodeLimit);

    static BDD
    bvec_gte(const Bvec& left, const Bvec& right, bool precise);

    static BDD
    bvec_gte_reduced(const Bvec& left, const Bvec& right, traverse_heuristic heu, unsigned int nodeLimit);

    static BDD
    bvec_slth(const Bvec& left, const Bvec& right, bool precise) {
        Cudd& manager = check_same_cudd(*left.m_manager, *right.m_manager);

        if (left.bitnum() == 0 || right.bitnum() == 0) {
            return manager.bddZero();
        }

        size_t size = left.bitnum() - 1;

        BDD differentSigns = precise ? left[size].AndP(~right[size]) : (left[size] & (~right[size]));
        if (differentSigns.IsOne())
        {
            // negative < positive
            return differentSigns;
        }
        else if (left[size].IsZero() && right[size].IsOne())
        {
            // positive < negative
            return manager.bddZero();
        }
        else
        {
            const Bvec &l_short = left.bvec_coerce(size);
            const Bvec &r_short = right.bvec_coerce(size);
            BDD equalSigns = precise ? left[size].XnorP(right[size]) : left[size].Xnor(right[size]);
            if (equalSigns.IsZero())    // don't need to compute lth which is possibly expensive
                return differentSigns;
            
            return precise
                ? differentSigns.OrP(equalSigns.AndP(bvec_lth(l_short, r_short, precise)))
                : differentSigns | (equalSigns & bvec_lth(l_short, r_short, precise));
        }
    }

    static BDD
    bvec_slth_reduced(const Bvec& left, const Bvec& right, traverse_heuristic heu, unsigned int nodeLimit) {
        Cudd& manager = check_same_cudd(*left.m_manager, *right.m_manager);

        if (left.bitnum() == 0 || right.bitnum() == 0) {
            return manager.bddZero();
        }

        size_t size = left.bitnum() - 1;

        BDD differentSigns = left[size].AndLim(~right[size], heu, nodeLimit);
        if (differentSigns.IsOne())
        {
            // negative < positive
            return differentSigns;
        }
        else if (left[size].IsZero() && right[size].IsOne())
        {
            // positive < negative
            return manager.bddZero();
        }
        else
        {
            const Bvec &l_short = left.bvec_coerce(size);
            const Bvec &r_short = right.bvec_coerce(size);
            BDD equalSigns = left[size].XnorLim(right[size], heu, nodeLimit);
            if (equalSigns.IsZero()) {    // don't need to compute lth which is possibly expensive
                return differentSigns;
            
            }
            return differentSigns.OrLim(equalSigns.AndLim(
                bvec_lth_reduced(l_short, r_short, heu, nodeLimit), heu, nodeLimit), heu, nodeLimit);
        }
    }

    static BDD
    bvec_slte(const Bvec& left, const Bvec& right, bool precise) {
        Cudd& manager = check_same_cudd(*left.m_manager, *right.m_manager);
        if (left.bitnum() == 0 || right.bitnum() == 0) {
            return manager.bddZero();
        }

        size_t size = left.bitnum() - 1;

        BDD differentSigns = precise ? left[size].AndP(~right[size]) : (left[size] & (~right[size]));
        if (differentSigns.IsOne())
        {
            // negative <= positive
            return differentSigns;
        }
        else if (left[size].IsZero() && right[size].IsOne())
        {
            // positive <= negative
            return manager.bddZero();
        }
        else
        {
            const Bvec &l_short = left.bvec_coerce(size);
            const Bvec &r_short = right.bvec_coerce(size);
            BDD equalSigns = precise ? left[size].XnorP(right[size]) : left[size].Xnor(right[size]);
            if (equalSigns.IsZero())    // don't need to compute lte which is possibly expensive
                return differentSigns;
            
            return precise
                ? differentSigns.OrP(equalSigns.AndP(bvec_lte(l_short, r_short, precise)))
                : differentSigns | (equalSigns & bvec_lte(l_short, r_short, precise));
        }
    }

    static BDD
    bvec_slte_reduced(const Bvec& left, const Bvec& right, traverse_heuristic heu, unsigned int nodeLimit) {
        Cudd& manager = check_same_cudd(*left.m_manager, *right.m_manager);
        
        if (left.bitnum() == 0 || right.bitnum() == 0) {
            return manager.bddZero();
        }

        size_t size = left.bitnum() - 1;

        BDD differentSigns = left[size].AndLim(~right[size], heu, nodeLimit);
        if (differentSigns.IsOne())
        {
            // negative <= positive
            return differentSigns;
        }
        else if (left[size].IsZero() && right[size].IsOne())
        {
            // positive <= negative
            return manager.bddZero();
        }
        else
        {
            const Bvec &l_short = left.bvec_coerce(size);
            const Bvec &r_short = right.bvec_coerce(size);
            BDD equalSigns = left[size].XnorLim(right[size], heu, nodeLimit);
            if (equalSigns.IsZero()) {    // don't need to compute lte which is possibly expensive
                return differentSigns;
            }
            
            return differentSigns.OrLim(equalSigns.AndLim(
                bvec_lte_reduced(l_short, r_short, heu, nodeLimit), heu, nodeLimit), heu, nodeLimit);
        }
    }

    static BDD
    bvec_sgth(const Bvec& left, const Bvec& right, bool precise);

    static BDD
    bvec_sgth_reduced(const Bvec& left, const Bvec& right, traverse_heuristic heu, unsigned int nodeLimit);

    static BDD
    bvec_sgte(const Bvec& left, const Bvec& right, bool precise);

    static BDD
    bvec_sgte_reduced(const Bvec& left, const Bvec& right, traverse_heuristic heu, unsigned int nodeLimit);

    static BDD
    bvec_equ(const Bvec& left, const Bvec& right, bool precise) {
       Cudd& manager = check_same_cudd(*left.m_manager, *right.m_manager);
       BDD p = manager.bddOne();

       if (left.bitnum() == 0 || right.bitnum() == 0 || left.bitnum() != right.bitnum()) {
           return manager.bddZero();
       }

       for (size_t i = 0U; i < left.bitnum(); ++i) {
           p = precise ? p.AndP(left[i].XnorP(right[i])) : p & left[i].Xnor(right[i]);
           if (p.IsZero())
           {
               return p;
           }
       }
       return p;
    }
    
    static BDD
    bvec_equ_reduced(const Bvec& left, const Bvec& right, traverse_heuristic heu, unsigned int nodeLimit) {
       Cudd& manager = check_same_cudd(*left.m_manager, *right.m_manager);
       BDD p = manager.bddOne();

       if (left.bitnum() == 0 || right.bitnum() == 0 || left.bitnum() != right.bitnum()) {
           return manager.bddZero();
       }

       for (size_t i = 0U; i < left.bitnum(); ++i) {
           p = p.AndLim(left[i].XnorLim(right[i], heu, nodeLimit), heu, nodeLimit);
           if (p.IsZero())
           {
               return p;
           }
       }
       return p;
    }

    static BDD
    bvec_nequ(const Bvec& left, const Bvec& right, bool precise) {
        Cudd& manager = check_same_cudd(*left.m_manager, *right.m_manager);
        BDD p = manager.bddZero();

        if (left.bitnum() == 0 || right.bitnum() == 0 || left.bitnum() != right.bitnum()) {
            return manager.bddZero();
        }

        for (size_t i = 0U; i < left.bitnum(); ++i) {
            p = precise ? p.OrP(left[i].XorP(right[i])) : p | left[i].Xor(right[i]);
            if (p.IsOne())
            {
                return p;
            }
        }
        return p;
    }

    static BDD
    bvec_nequ_reduced(const Bvec& left, const Bvec& right, traverse_heuristic heu, unsigned int nodeLimit) {
        Cudd& manager = check_same_cudd(*left.m_manager, *right.m_manager);
        BDD p = manager.bddZero();

        if (left.bitnum() == 0 || right.bitnum() == 0 || left.bitnum() != right.bitnum()) {
            return manager.bddZero();
        }

        for (size_t i = 0U; i < left.bitnum(); ++i) {
            p = p.OrLim(left[i].XorLim(right[i], heu, nodeLimit), heu, nodeLimit);
            if (p.IsOne())
            {
                return p;
            }
        }
        return p;
    }

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
    operator<<(int con) const { return bvec_shlfixed(con, m_manager->bddZero()); }

    Bvec
    operator<<(const Bvec& other) const { return bvec_shl(*this, other, m_manager->bddZero(), false); }

    Bvec
    operator>>(int con) const { return bvec_shrfixed(con, m_manager->bddZero()); }

    Bvec
    operator>>(const Bvec& other) const { return bvec_shr(*this, other, m_manager->bddZero(), false); }

    Bvec
    operator+(const Bvec& other) const { return bvec_add(*this, other, false); }

    Bvec
    operator+=(const Bvec& other) { *this = bvec_add(*this, other, false); return *this; }


    Bvec
    operator-(const Bvec& other) { return bvec_sub(*this, other, false); }

    Bvec
    operator-=(const Bvec& other) { *this = bvec_sub(*this, other, false); return *this; }

    Bvec
    operator*(int con) const { return bvec_mulfixed(con, false); }

    Bvec
    operator*=(int con) { this->bvec_mulfixed(con, false); return *this; }

    Bvec
    operator*(const Bvec& other) const { return bvec_mul(*this, other, false); }

    Bvec
    operator*=(const Bvec& other) { *this = bvec_mul(*this, other, false); return *this; }

    BDD
    operator<(const Bvec& other) const { return bvec_lth(*this, other, false); }

    BDD
    operator<=(const Bvec& other) const { return bvec_lte(*this, other, false); }

    BDD
    operator>(const Bvec& other) const { return bvec_gth(*this, other, false); }

    BDD
    operator>=(const Bvec& other) const { return bvec_gte(*this, other, false); }

    BDD
    operator==(const Bvec& other) const { return bvec_equ(*this, other, false); }

    BDD
    operator!=(const Bvec& other) const { return !(*this == other); }

    unsigned int bddNodes()
    {
        auto count = 0U;

        for (const auto &bdd : m_bitvec)
        {
            count += bdd.nodeCount();
        }

        return count;
    }

    bool isPrecise() const
    {
        for (const auto &bdd : m_bitvec)
        {
            if (bdd.IsUnknown())
            {
                return false;
            }
        }

        return true;
    }

private:

    static Cudd&
    check_same_cudd(Cudd& first, Cudd& second);

    static void
    bvec_div_rec(Bvec& divisor, Bvec& remainder, Bvec& result, size_t step, bool precise);

    static BDD
    bdd_and(const BDD& first, const BDD& second);

    static BDD
    bdd_xor(const BDD& first, const BDD& second);

    static BDD
    bdd_or(const BDD& first, const BDD& second);

    static BDD
    bdd_not(const BDD& src);

    void
    swap(Bvec& other);

    static Bvec
    reserve(Cudd& manager, size_t bitnum);

    static void
    reserve(Bvec& bitvector, size_t bitnum);
};

} // cudd

#endif //BDD_BVEC_H
