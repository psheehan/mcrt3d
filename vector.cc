#ifndef VECTOR_CC
#define VECTOR_CC

#include <stdio.h>
#include <cmath>

template<typename Type, int D>
struct Vector {
    Type v[D];

    Vector() {
        for (int i=0; i<D; i++) v[i] = 0;
    }

    Vector(Type scalar) {
        for (int i=0; i<D; i++) v[i] = scalar;
    }

    Vector(Type x, Type y) {
        v[0] = x;
        v[1] = y;
    }

    Vector(Type x, Type y, Type z) {
        v[0] = x;
        v[1] = y;
        v[2] = z;
    }

    Type operator[] ( int i ) const {
        return v[i];
    }

    Type& operator[] ( int i ) {
        return *(v+i);
    }

    Vector<Type, D> operator+() const {
        return *this;
    }

    Vector<Type, D> operator-() const {
        return Vector<Type, D>(-v[0],-v[1],-v[2]);
    }

    template<typename Type2>
    Vector<Type, D> operator= (Type2 scalar) {
        for (int i=0; i<D; i++) v[i] = scalar;
        return *this;
    }

    template<typename Type2>
    Vector<Type, D> operator= (Vector<Type2, D> rhs) {
        for (int i=0; i<D; i++) v[i] = rhs[i];
        return *this;
    }

    template<typename Type2>
    Vector<Type, D> operator+= (Type2 scalar) {
        for (int i=0; i<D; i++) v[i] += scalar;
        return *this;
    }

    template<typename Type2>
    Vector<Type, D> operator+= (const Vector<Type2, D> rhs) {
        for (int i=0; i<D; i++) v[i] += rhs[i];
        return *this;
    }

    template<typename Type2>
    Vector<Type, D> operator-= (Type2 scalar) {
        for (int i=0; i<D; i++) v[i] -= scalar;
        return *this;
    }

    template<typename Type2>
    Vector<Type, D> operator-= (const Vector<Type2, D> rhs) {
        for (int i=0; i<D; i++) v[i] -= rhs[i];
        return *this;
    }

    template<typename Type2>
    Vector<Type, D> operator*= (Type2 scalar) {
        for (int i=0; i<D; i++) v[i] *= scalar;
        return *this;
    }

    template<typename Type2>
    Vector<Type, D> operator*= (const Vector<Type2, D> rhs) {
        for (int i=0; i<D; i++) v[i] *= rhs[i];
        return *this;
    }

    template<typename Type2>
    Vector<Type, D> operator/= (Type2 scalar) {
        for (int i=0; i<D; i++) v[i] /= scalar;
        return *this;
    }

    template<typename Type2>
    Vector<Type, D> operator/= (const Vector<Type2, D> rhs) {
        for (int i=0; i<D; i++) v[i] /= rhs[i];
        return *this;
    }

    Type norm() const {
        Type result = 0;
        for (int i=0; i<D; i++) result += v[i]*v[i];
        return sqrt(result);
    }

};

template<typename Type1, typename Type2, int D>
Vector<Type1, D> operator+ (const Vector<Type1, D> lhs, const Vector<Type2, D> 
        rhs) {
    Vector<Type1, D> new_vector(0);
    for (int i=0; i<D; i++) new_vector.v[i] = lhs.v[i] + rhs.v[i];
    return new_vector;
}

template<typename Type1, typename Type2, int D>
Vector<Type1, D> operator+ (const Vector<Type1, D> lhs, Type2 rhs) {
    Vector<Type1, D> new_vector(0);
    for (int i=0; i<D; i++) new_vector.v[i] = lhs.v[i] + rhs;
    return new_vector;
}

template<typename Type1, typename Type2, int D>
Vector<Type1, D> operator+ (Type1 lhs, const Vector<Type2, D> rhs) {
    Vector<Type1, D> new_vector(0);
    for (int i=0; i<D; i++) new_vector.v[i] = lhs + rhs.v[i];
    return new_vector;
}

template<typename Type1, typename Type2, int D>
Vector<Type1, D> operator- (const Vector<Type1, D> lhs, const Vector<Type2, D> 
        rhs) {
    Vector<Type1, D> new_vector(0);
    for (int i=0; i<D; i++) new_vector.v[i] = lhs.v[i] - rhs.v[i];
    return new_vector;
}

template<typename Type1, typename Type2, int D>
Vector<Type1, D> operator- (const Vector<Type1, D> lhs, Type2 rhs) {
    Vector<Type1, D> new_vector(0);
    for (int i=0; i<D; i++) new_vector.v[i] = lhs.v[i] - rhs;
    return new_vector;
}

template<typename Type1, typename Type2, int D>
Vector<Type1, D> operator- (Type1 lhs, const Vector<Type2, D> rhs) {
    Vector<Type1, D> new_vector(0);
    for (int i=0; i<D; i++) new_vector.v[i] = lhs - rhs.v[i];
    return new_vector;
}

template<typename Type1, typename Type2, int D>
Type1 operator* (const Vector<Type1, D> lhs, const Vector<Type2, D> rhs) {
    Type1 result = 0;
    for (int i=0; i<D; i++) result += lhs.v[i] * rhs.v[i];
    return result;
}

template<typename Type1, typename Type2, int D>
Vector<Type1, D> operator* (const Vector<Type1, D> lhs, Type2 rhs) {
    Vector<Type1, D> new_vector(0);
    for (int i=0; i<D; i++) new_vector.v[i] = lhs.v[i] * rhs;
    return new_vector;
}

template<typename Type1, typename Type2, int D>
Vector<Type1, D> operator* (Type1 lhs, const Vector<Type2, D> rhs) {
    Vector<Type1, D> new_vector(0);
    for (int i=0; i<D; i++) new_vector.v[i] = lhs * rhs.v[i];
    return new_vector;
}

template<typename Type1, typename Type2, int D>
Vector<Type1, D> operator/ (const Vector<Type1, D> lhs, Type2 rhs) {
    Vector<Type1, D> new_vector(0);
    for (int i=0; i<D; i++) new_vector.v[i] = lhs.v[i] / rhs;
    return new_vector;
}

template<typename Type1, typename Type2, int D>
Vector<Type1, D> operator/ (Type1 lhs, const Vector<Type2, D> rhs) {
    Vector<Type1, D> new_vector(0);
    for (int i=0; i<D; i++) new_vector.v[i] = lhs / rhs.v[i];
    return new_vector;
}

#endif
