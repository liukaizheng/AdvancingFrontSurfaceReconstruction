#ifndef PLAUSIBILITY_H
#define PLAUSIBILITY_H

#include <set>

template <typename Scalar, typename Index>
struct Plausibility {
    Plausibility(const Scalar& _p, const Index& _f, const Index& _ei):p(_p), f(_f), ei(_ei){}
    Scalar p;
    Index f;
	Index ei;
};

template <typename Scalar, typename Index>
struct std::less<Plausibility<Scalar, Index>*> {
    typedef Plausibility<Scalar, Index> Base;
    bool operator()(const Base* p1, const Base* p2) const
    {
        if(p1->p < p2->p)
            return true;
        else if(p1->p > p2->p)
            return false;
        else if(p1->f < p2->f)
            return true;
        else if(p1->f > p2->f)
            return false;
        else if(p1->ei < p2->ei)
            return true;
        else if(p1->ei > p2->ei)
            return false;
        else
            return false;
    }
};

#endif
