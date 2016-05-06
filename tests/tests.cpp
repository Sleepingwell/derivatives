#include <stdio.h>
#include "../include/SKT/derivatives/derivatives.hpp"

#define N(n) "f"#n" = %g\n\n"

// A and B are structs that hold a single member 'variable'
struct A {
	SK_DERIV_MEMBER_VARIABLE(a);
	A(double a) : SK_DERIV_INIT_VARIABLE(a, a) {}
};

struct B {
	SK_DERIV_MEMBER_VARIABLE(b);
	B(double b) : SK_DERIV_INIT_VARIABLE(b, b) {}
};

// C is a struct that holds a few member 'variables'
struct C {
	SK_DERIV_MEMBER_VARIABLE(a);
	SK_DERIV_MEMBER_VARIABLE(b);
	SK_DERIV_MEMBER_VARIABLE(c);
	C(void) : SK_DERIV_INIT_VARIABLE(a, 2.0), SK_DERIV_INIT_VARIABLE(b, 3.0), SK_DERIV_INIT_VARIABLE(c, 4.0) {}

	SK_DERIV_DECLARE_DERIVATIVES(C, (2.0 * pow(a, b*c) + a*b*c*c*c / 5.4))
};

// F and G are structs that hold a few member variables and a (differentiable)
// class or member 'variable' of a class (since the implementation is the same, G is redundant
// as classes produced for either template are specific to the type EXT... but you knew this).
template<typename EXT>
struct F {
	SK_DERIV_MEMBER_VARIABLE(a);
	SK_DERIV_MEMBER_VARIABLE(b);
	SK_DERIV_MEMBER_VARIABLE(c);
	F(double a, double b, double c, typename boost::call_traits<EXT>::param_type other) :
	// note that normal declarations are also possible (i.e. that SK_DERIV_INIT_VARIABLE is just for convenience).
		a(a), b(b), c(c), other(other) {}

	SK_DERIV_DECLARE_DERIVATIVES(F, (other * pow(a, b*c) + a*b*pow(c, 3) / 3.9))

private:
	typename boost::call_traits<EXT>::const_reference other;
};

template<typename EXT>
struct G {
	SK_DERIV_MEMBER_VARIABLE(a);
	SK_DERIV_MEMBER_VARIABLE(b);
	SK_DERIV_MEMBER_VARIABLE(c);
	G(double a, double b, double c, typename boost::call_traits<EXT>::param_type other) :
		SK_DERIV_INIT_VARIABLE(a, a), SK_DERIV_INIT_VARIABLE(b, b), SK_DERIV_INIT_VARIABLE(c, c), other(other) {}

	//SK_DERIV_DECLARE_DERIVATIVES(G, (other.c * pow(a, b*c) + a*b*pow(c, 3) / 3.9))
	SK_DERIV_DECLARE_DERIVATIVES(G, (other))

private:
	typename boost::call_traits<EXT>::const_reference other;
};

// a few Global 'variables'.
SK_DERIV_GLOBAL_VARIABLE(a, 10.0);
SK_DERIV_GLOBAL_VARIABLE(b, 11.0);
SK_DERIV_GLOBAL_VARIABLE(c, 12.0);
SK_DERIV_GLOBAL_VARIABLE(d, 13.0);

#define MODEL 5.0*a + b / c
// P(a) for 'print a' (symbolic)
// D(a, b) da/db
// E(a) for 'evaluate a'
int main(int argc, char* argv[]){
	//double[] vect = {
	//	E(D(D(MODEL, a), b)),
	//	E(D(MODEL, b))
	//}
//
//	printf("\
//-----------------------------------------------------------------\n\
//1st to 4th order derivatives of an arbitrary expression.\n \
//-----------------------------------------------------------------\n");
	printf("1st order\n");
	P(D((c*a + b) / ((5.0*a + b / c) * d), c));
	printf("2nd order\n");
	P(D(D((c*a + b) / ((5.0*a + b / c) * d), c), c));
	//printf("3rd order\n");
	//P(D(D(D((c*a + b) / ((5.0*a + b / c) * d), c), c), c));
	//printf("4th order\n");
	//P(D(D(D(D((c*a + b) / ((5.0*a + b / c) * d), c), c), c), c));
    //printf("5th order\n");
	//P(D(D(D(D(D((c*a + b) / ((5.0*a + b / c) * d), c), c), c), c), c)); //this blows VC++ heap, but worked on gcc.
//
	printf("\
-----------------------------------------------------------------\n\
1st to 4th order derivatives of another arbitrary expression.\n \
-----------------------------------------------------------------\n");
	printf("1st order\n");
	P(D(a/c, c));
	//printf("2nd order\n");
	//P(D(D(a/c, c), c));
	//printf("3rd order\n");
	//P(D(D(D(a/c, c), c), c));
//	printf("4th order\n");
	//P(D(D(D(D(a/c, c), c), c), c));
//	////P(D(D(D(D(D((c*a + b) / ((5.0*a + b / c) * d), c), c), c), c), c)); //this blows VC++ heap, but worked on gcc.
//
	printf("\
-----------------------------------------------------------------\n \
A few misc derivatives (easy to see what simplification does).\n \
-----------------------------------------------------------------\n");
	P(D(a*b*b, b) + D(c*b*b, b));

	P(D(exp(b), b));

	P(D(pow(EXP, b), b));

	printf("\n\n\
//-----------------------------------------------------------------\n \
//Derivitives of 'simple' classes.\n \
//-----------------------------------------------------------------\n");
//	C cc;
//	E(cc);
//	printf("d11 = %g\n", E(c*D(cc, cc.a)));
//	printf("d12 = %g\n", E(c*D(D(cc, cc.a), cc.b)));
//	printf("d13 = %g\n", E(c*D(D(D(cc, cc.a), cc.b), cc.c)));
////	//printf("d14 = %g\n", E(c*D(D(D(D(cc, cc.a), cc.b), cc.c), cc.c)));
//
//	printf("\n");
//
//	A Aa(42.0);
//	B Bb(7.1);
//	P(D(D(pow(Aa.a, 2.0) + exp(Bb.b), Aa.a), Bb.b));
//	P(D(D(pow(Aa.a, 2.0) * exp(Bb.b), Aa.a), Bb.b));
//
//	printf("\n\n\
//-----------------------------------------------------------------\n \
//Derivitives of classes that hold references to variables from other classes.\n \
//-----------------------------------------------------------------\n");
//	F<C::typeofc> fc(1.0, 2.0, 3.0, cc.c);
//	printf("d21 = %g\n", E(fc));
//	printf("d22 = %g\n", E(c*D(fc, cc.c)));
//	printf("d23 = %g\n", E(fc*D(fc + cc, cc.c)));
//	//printf("d24 = %g\n", E(D(pow(fc, 4), cc.c)));
//
//	printf("\n\n\
//-----------------------------------------------------------------\n \
//Derivitives of classes that hold references to other classes.\n \
//-----------------------------------------------------------------\n");
//	G<C> gc(1.0, 2.0, 3.0, cc);
//	printf("d21 = %g\n", E(gc));
//	printf("d22 = %g\n", E(c*D(gc, cc.c)));
//	printf("d23 = %g\n", E(gc*D(gc + cc, cc.c)));
//	//printf("d24 = %g\n", E(D(pow(gc, 4), cc.c)));
//
////	printf("\n\n\
////-----------------------------------------------------------------\n \
////test the chain rule \
////-----------------------------------------------------------------\n");
////	printf("%g == %g\n", E(D(gc, cc.c)), 0.0); //E(D(gc, cc) * D(cc, cc.c)));
////	printf("%g\n", E(D(gc, cc))); //E(D(gc, cc) * D(cc, cc.c)));

	return 0;
}
