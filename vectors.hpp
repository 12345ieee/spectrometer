#ifndef _VECTORS_HPP
#define _VECTORS_HPP

#include <iostream>

#include "TVector2.h"
#include "TVector3.h"

std::ostream& operator<<(std::ostream& ost, TVector2& vec)
{
    return ost << vec.X() << " " << vec.Y();
}

std::ostream& operator<<(std::ostream& ost, TVector3& vec)
{
    return ost << vec.X() << " " << vec.Y() << " " << vec.Z();
}

std::istream& operator>>(std::istream& ist, TVector2& vec)
{
    double x, y;
    ist >> x >> y;
    vec = TVector2(x, y);
    return ist;
}

std::istream& operator>>(std::istream& ist, TVector3& vec)
{
    double x, y, z;
    ist >> x >> y >> z;
    vec = TVector3(x, y, z);
    return ist;
}

#endif // _VECTORS
