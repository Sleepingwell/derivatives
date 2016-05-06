/*
 * file derivatives.hpp
 * Copyright (C) 2008-2011       Simon Knapp
 *
 * This progrm is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SK_DERIVIITIVES_HEADER_INCLUDED_298YWDOIKJHWOTUYHWOKJGNLIUHWT
#define SK_DERIVIITIVES_HEADER_INCLUDED_298YWDOIKJHWOTUYHWOKJGNLIUHWT

#include <cmath>
#include <boost/static_assert.hpp>
#include <boost/call_traits.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/or.hpp>
#include <boost/type_traits/add_const.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/remove_cv.hpp>
#include <boost/type_traits/remove_reference.hpp>

#define EXP 2.7182818284590452354
#define DTC(D_PARAM_T) typename SKT::derivatives::type_chooser< D_PARAM_T >::bare_type
#define DPARM(D_PARAM_T) typename boost::call_traits< D_PARAM_T >::param_type
#define DPARMET(D_PARAM_T) DPARM(D_PARAM_T)
#define DGD(D_STRUCT_T, D_STRUCT_V, D_WRT_T, D_WRT_V ) D_STRUCT_T::template GetDeravitive<D_WRT_T>::deriv(D_STRUCT_V, D_WRT_V)

namespace SKT { 
	namespace derivatives { // forward defs of types needed by operators.
		template<typename L, typename OP, typename R> class Creator;
		template<typename IN> struct type_chooser;
		struct Numeric;
		struct OpPlus; struct OpMinus; struct OpTimes; struct OpDivide; struct OpPower; struct OpLog;
	}
}
//--------------------------------------------------------------------
// Operator forward defs
//--------------------------------------------------------------------
template<typename L, typename R>
typename SKT::derivatives::Creator<L, SKT::derivatives::OpPlus, R>::type operator+(L const& l, R const& r);

template<typename L, typename R>
typename SKT::derivatives::Creator<L, SKT::derivatives::OpMinus, R>::type operator-(L const& l, R const& r);

template<typename L, typename R>
typename SKT::derivatives::Creator<L, SKT::derivatives::OpTimes, R>::type operator*(L const& l, R const& r);

template<typename L, typename R>
typename SKT::derivatives::Creator<L, SKT::derivatives::OpDivide, R>::type operator/(L const& l, R const& r);

template<typename L, typename R>
typename SKT::derivatives::Creator<L, SKT::derivatives::OpPower, R>::type pow(L const& l, R const& r);

template<typename R>
typename SKT::derivatives::Creator<SKT::derivatives::Numeric, SKT::derivatives::OpPower, R>::type exp(R const& r);

template<typename L>
typename SKT::derivatives::Creator<L, SKT::derivatives::OpLog, SKT::derivatives::Numeric>::type log(L const& l);


namespace SKT { namespace derivatives {

	using namespace boost;

	struct differentiable_class_tag {};

	struct Numeric; // forward def
	template<typename L, typename OP, typename R> struct Expression; // forward def
	template<int N, typename T, typename WRT> struct DerivativeProxy; // forward def

	struct Zero : differentiable_class_tag {
		static double eval(void) const { return 0.0; }
		template<typename T> struct GetDeravitive { static Zero deriv(DPARM(Zero), DPARM(T)) { return Zero(); } };
		static void print(void) const { printf("Zero"); }
	};

	struct One : differentiable_class_tag {
		static double eval(void) const { return 1.0; }
		template<typename T> struct GetDeravitive { static Zero deriv(DPARM(One), DPARM(T)) { return Zero(); } };
		static void print(void) const { printf("One"); }
	};

	struct Numeric {
		Numeric(double t) : val_(t) {}
		Numeric(Numeric const& rhs) : val_(rhs.val_) {}
		double eval(void) const { return val_; }
		template<typename T> struct GetDeravitive { static Zero deriv(DPARM(Numeric), DPARM(T)) { return Zero(); } };
		void print(void) const { printf("%g", val_); }

	private:
		double val_;
	};

//--------------------------------------------------------------------
// choose types to be held in an expression
//--------------------------------------------------------------------
	template< typename T >
	struct remove_all : remove_cv< typename remove_reference< T >::type > {};

	template<typename IN>
	struct type_chooser {
		typedef typename remove_all< IN >::type raw_type;
		typedef typename mpl::if_< is_arithmetic< raw_type >, Numeric, raw_type >::type bare_type;

		typedef typename mpl::if_<
			is_same< Numeric, bare_type >,
			typename add_const< Numeric >::type, // keep arithmetic types and Numerics as const Numeric values.
			typename mpl::if_<
				is_base_and_derived< differentiable_class_tag, bare_type >,
				typename add_const< bare_type >::type, // keep internal types as values.
				typename call_traits< bare_type >::const_reference // keep everything else as const references.
			>::type
		>::type held_type;
	};

//--------------------------------------------------------------------
// get return types
//--------------------------------------------------------------------
	template<typename LB, typename OP, typename RB, typename L, typename R>
	struct CreatorImpl : OP::Simplify::template apply<LB, RB, L, R> {}; //struct CreatorImpl : Simplify< mpl::apply< typename OP::Simplify, LB, RB, L, R> > {}; //slows compile down significantly

	template<typename L, typename OP, typename R>
	struct CreatorImpl<Numeric, OP, Numeric, L, R> {
		typedef Numeric type;
		static type create(DPARM(L) l, DPARM(R) r) {
			return Numeric(OP::eval<L, R>(l.eval(), r.eval()));
		}
	};

	//template<typename OP, typename L, typename R, typename T>
	//struct CreatorImpl<T, OP, double, L, R> : CreatorImpl<T, OP, Numeric, L, Numeric> {
	//	typedef CreatorImpl<T, OP, Numeric, L, Numeric> Base;
	//	static typename Base::type create(DPARM(L) l, DPARM(R) r) {
	//		return Base::create(l, Numeric(r));
	//	}
	//};

	//template<typename OP, typename L, typename R, typename T>
	//struct CreatorImpl<double, OP, T, L, R> : CreatorImpl<Numeric, OP, T, Numeric, R> {
	//	typedef CreatorImpl<Numeric, OP, T, Numeric, R> Base;
	//	static typename Base::type create(DPARM(L) l, DPARM(R) r) {
	//		return Base::create(Numeric(l), r);
	//	}
	//};

	//template<typename OP, typename L, typename R>
	//struct CreatorImpl<double, OP, double, L, R> {
	//	typedef Numeric type;
	//	static type create(double l, double r) {
	//		return Numeric(OP::eval(l, r));
	//	}
	//};

	template<typename L, typename OP, typename R>
	class Creator {
		typedef DTC(L) left_bare_type;
		typedef DTC(R) right_bare_type;
	public:
		typedef typename CreatorImpl< left_bare_type, OP, right_bare_type, L, R >::type type;
		static type create(DPARM(L) l, DPARM(R) r) {
			return CreatorImpl< left_bare_type, OP, right_bare_type, L, R >::create(l, r);
		}
	};

//--------------------------------------------------------------------
// derivative types
//--------------------------------------------------------------------
	template<typename T, typename WRT>
	struct DerivatorImpl {
		typedef typename T::template DerivativeType< WRT >::type type;
	};

	template<typename T>
	struct DerivatorImpl<T, T> {
		typedef One type;
	};

	//template<typename T>
	//struct DerivatorImpl<double, T> {
	//	typedef Zero type;
	//};

	template<typename T>
	struct DerivatorImpl<Numeric, T> {
		typedef Zero type;
	};

	template<typename T>
	struct DerivatorImpl<Zero, T> {
		typedef Zero type;
	};

	template<typename T>
	struct DerivatorImpl<One, T> {
		typedef Zero type;
	};

	template<int N, typename T, typename WRTP, typename WRT>
	struct DerivatorImpl< DerivativeProxy<N, T, WRTP>, WRT > {
		typedef DerivativeProxy<N+1, DerivativeProxy<N, T, WRTP>, WRT> type;
	};

	template<typename L, typename OP, typename R, typename WRT>
	struct DerivatorImpl< Expression<L, OP, R>, WRT > {
		typedef typename OP::DerivType::template apply< typename Expression<L, OP, R>::left_type, typename Expression<L, OP, R>::right_type, WRT >::type type;
	};

	template<typename T, typename WRT>
	struct Derivator : DerivatorImpl< DTC(T), DTC(WRT) > {
	};

//--------------------------------------------------------------------
// utilities to get the types for variables.
//--------------------------------------------------------------------
	template<typename T, typename WRT >
	struct DerivativeType { 
		typedef Zero type; 
	};

	template<typename T>
	struct DerivativeType< T, T > {
		typedef SKT::derivatives::One type;
	};

//--------------------------------------------------------------------
// core type
//--------------------------------------------------------------------
	template<typename L, typename OP, typename R>
	struct Expression : differentiable_class_tag {
	private:
		typedef type_chooser< L > left_tc_type;
		typedef type_chooser< R > right_tc_type;
		typedef typename left_tc_type::held_type held_left_type;
		typedef typename right_tc_type::held_type held_right_type;
		
	public:
		typedef typename left_tc_type::bare_type left_type;
		typedef typename right_tc_type::bare_type right_type;

		Expression(DPARM(L) l, DPARM(R) r) : l(l), r(r) {
		}

		Expression(DPARM(Expression) rhs) : l(rhs.l), r(rhs.r) {
		}

		double eval(void) const {
			return OP::eval<held_left_type, held_right_type>(l, r);
		}

		template<typename WRT>
		struct GetDeravitive {
			static typename Derivator<Expression, WRT>::type deriv(DPARM(Expression) expr, DPARM(WRT) wrt) {
				return OP::deriv<WRT, left_type, right_type>(wrt, expr.l, expr.r);
			};
		};

		void print(void) const {
			OP::print<left_type, right_type>(l, r);
		}

		held_left_type l;
		held_right_type r;
	};

//--------------------------------------------------------------------
// operators
//--------------------------------------------------------------------
	namespace Private {
		template<typename L, typename OP, typename R, typename LB, typename RB>
		struct DefaultExpression {
			typedef Expression<LB, OP, RB> type;
			static type create(DPARM(L) l, DPARM(R) r) {
				return type(l, r);
			}
		};
	} // end namespace Private

	#define SK_DERIV_SIMPLIFY_BASE(OP) template<typename LB, typename RB, typename L, typename R>  struct apply : Private::DefaultExpression<L, OP, R, LB, RB> {}

// +
	struct OpPlus {
		template<typename L, typename R>
		static double eval(DPARM(L) left, DPARM(R) right) {
			return left.eval() + right.eval();
		}

		struct DerivType {
			template<typename L, typename R, typename WRT>
			struct apply : Creator< typename Derivator<L, WRT>::type, OpPlus, typename Derivator<R, WRT>::type > {};
		};

		struct Simplify {
			SK_DERIV_SIMPLIFY_BASE(OpPlus);

			template<typename RB, typename L, typename R>
			struct apply<Zero, RB, L, R> {
				typedef typename type_chooser< R >::held_type type;
				static type create(DPARM(L) l, DPARM(R) r) { return r; }
			};

			template<typename LB, typename L, typename R>
			struct apply<LB, Zero, L, R> {
				typedef typename type_chooser< L >::held_type type;
				static type create(DPARM(L) l, DPARM(R) r) { return l; }
			};

			template<typename L, typename R>
			struct apply<Zero, Zero, L, R> {
				typedef Zero type;
				static type create(DPARM(L) l, DPARM(R) r) { return Zero(); }
			};
		};

		template<typename WRT, typename L, typename R>
		static typename DerivType::template apply<L, R, WRT>::type deriv(DPARM(WRT) wrt, DPARM(L) left, DPARM(R) right) {
			return DGD(L, left, WRT, wrt) + DGD(R, right, WRT, wrt);
		}

		template<typename L, typename R>
		static void print(DPARM(L) l, DPARM(R) r) {
			printf("("); l.print(); printf(" + "); r.print(); printf(")");
		}
	};

// -
	struct OpMinus {
		template<typename L, typename R>
		static double eval(DPARM(L) left, DPARM(R) right) {
			return left.eval() - right.eval();
		}

		struct DerivType {
			template<typename L, typename R, typename WRT>
			struct apply : Creator<
				typename Derivator<L, WRT>::type,
				OpMinus,
				typename Derivator<R, WRT>::type
			> {};
		};

		struct Simplify {
			SK_DERIV_SIMPLIFY_BASE(OpMinus);

			template<typename LB, typename L, typename R>
			struct apply<LB, Zero, L, R> {
				typedef typename type_chooser< L >::held_type type;
				static type create(DPARM(L) l, DPARM(R) r) { return l; }
			};

			template<typename L, typename R>
			struct apply<Zero, Zero, L, R> {
				typedef Zero type;
				static type create(DPARM(L) l, DPARM(R) r) { return Zero(); }
			};
		};

		template<typename WRT, typename L, typename R>
		static typename DerivType::template apply<L, R, WRT>::type deriv(WRT const& wrt, DPARM(L) left, DPARM(R) right) {
			return DGD(L, left, WRT, wrt) - DGD(R, right, WRT, wrt);
		}

		template<typename L, typename R>
		static void print(L const& l, R const& r) {
			printf("("); l.print(); printf(" - "); r.print(); printf(")");
		}
	};

// *
	struct OpTimes {
		template<typename L, typename R>
		static double eval(DPARM(L) left, DPARM(R) right) {
			return left.eval() * right.eval();
		}

		struct DerivType {
			template<typename L, typename R, typename WRT>
			struct apply : Creator<
				typename Creator< R, OpTimes, typename Derivator<L, WRT>::type >::type,
				OpPlus,
				typename Creator< L, OpTimes, typename Derivator<R, WRT>::type >::type
			> {};
		};

		struct Simplify {
			SK_DERIV_SIMPLIFY_BASE(OpTimes);

			template<typename RB, typename L, typename R>
			struct apply<Zero, RB, L, R> {
				typedef Zero type;
				static type create(DPARM(L) l, DPARM(R) r) { return Zero(); }
			};

			template<typename LB, typename L, typename R>
			struct apply<LB, Zero, L, R> {
				typedef Zero type;
				static type create(DPARM(L) l, DPARM(R) r) { return Zero(); }
			};

			template<typename RB, typename L, typename R>
			struct apply<One, RB, L, R> {
				typedef typename type_chooser< R >::held_type type;
				static type create(DPARM(L) l, DPARM(R) r) { return r; }
			};

			template<typename LB, typename L, typename R>
			struct apply<LB, One, L, R> {
				typedef typename type_chooser< L >::held_type type;
				static type create(DPARM(L) l, DPARM(R) r) { return l; }
			};

			template<typename L, typename R>
			struct apply<One, Zero, L, R> {
				typedef Zero type;
				static type create(DPARM(L) l, DPARM(R) r) { return Zero(); }
			};

			template<typename L, typename R>
			struct apply<Zero, One, L, R> {
				typedef Zero type;
				static type create(DPARM(L) l, DPARM(R) r) { return Zero(); }
			};

			template<typename L, typename R>
			struct apply<One, One, L, R> {
				typedef One type;
				static type create(DPARM(L) l, DPARM(R) r) { return One(); }
			};

			template<typename L, typename R>
			struct apply<Zero, Zero, L, R> {
				typedef Zero type;
				static type create(DPARM(L) l, DPARM(R) r) { return Zero(); }
			};
		};

		template<typename WRT, typename L, typename R>
		static typename DerivType::template apply<L, R, WRT>::type deriv(WRT const& wrt, DPARM(L) left, DPARM(R) right) {
			return right * DGD(L, left, WRT, wrt) + left * DGD(R, right, WRT, wrt);
		}

		template<typename L, typename R>
		static void print(L const& l, R const& r) {
			printf("("); l.print(); printf(" * "); r.print(); printf(")");
		}
	};

	struct OpDivide; // forward def

// log
	struct OpLog {
		template<typename L, typename R>
		static double eval(DPARM(L) left, DPARM(R) right) {
			return log(left.eval());
		}

		struct DerivType {
			template<typename L, typename R, typename WRT>
			struct apply : Creator< typename Derivator<L, WRT>::type, OpDivide, L > {};
		};

		struct Simplify {
			SK_DERIV_SIMPLIFY_BASE(OpLog);

			template<typename RB, typename L, typename R>
			struct apply<One, RB, L, R> {
				typedef Zero type;
				static type create(DPARM(L) l, DPARM(R) r) { return Zero(); }
			};
		};

		template<typename WRT, typename L, typename R>
		static typename DerivType::template apply<L, R, WRT>::type deriv(WRT const& wrt, DPARM(L) left, DPARM(R) right) {
			return DGD(L, left, WRT, wrt) / left;
		}

		template<typename L, typename R>
		static void print(L const& l, R const& r) {
			printf("log("); l.print(); printf(")");
		}
	};

// ^
	struct OpPower {
		template<typename L, typename R>
		static double eval(DPARM(L) left, DPARM(R) right) {
			return pow(left.eval(), right.eval());
		}

		struct DerivType {
			template<typename L, typename R, typename WRT>
			struct apply : Creator<
				typename Creator<
					typename Creator< typename Creator<L, OpLog, double>::type, OpTimes, typename Derivator<R, WRT>::type >::type,
					OpPlus,
					typename Creator< typename Derivator< typename Creator<L, OpLog, double>::type, WRT >::type, OpTimes, R >::type
				>::type,
				OpTimes,
				typename Creator< L, OpPower, R >::type
			> {};
		};

		struct Simplify {
			SK_DERIV_SIMPLIFY_BASE(OpPower);

			template<typename RB, typename L, typename R>
			struct apply<Zero, RB, L, R> {
				typedef Zero type;
				static type create(DPARM(L) l, DPARM(R) r) { return Zero(); }
			};

			template<typename RB, typename L, typename R>
			struct apply<One, RB, L, R> {
				typedef One type;
				static type create(DPARM(L) l, DPARM(R) r) { return One(); }
			};

			template<typename LB, typename L, typename R>
			struct apply<LB, Zero, L, R> {
				typedef One type;
				static type create(DPARM(L) l, DPARM(R) r) { return One(); }
			};

			template<typename LB, typename L, typename R>
			struct apply<LB, One, L, R> {
				typedef typename type_chooser< L >::held_type type;
				static type create(DPARM(L) l, DPARM(R) r) { return l; }
			};

			template<typename L, typename R>
			struct apply<Zero, Zero, L, R> {
				typedef One type;
				static type create(DPARM(L) l, DPARM(R) r) { return One(); }
			};

			template<typename L, typename R>
			struct apply<One, One, L, R> {
				typedef One type;
				static type create(DPARM(L) l, DPARM(R) r) { return One(); }
			};

			template<typename L, typename R>
			struct apply<Zero, One, L, R> {
				typedef Zero type;
				static type create(DPARM(L) l, DPARM(R) r) { return Zero(); }
			};

			template<typename L, typename R>
			struct apply<One, Zero, L, R> {
				typedef One type;
				static type create(DPARM(L) l, DPARM(R) r) { return One(); }
			};

			template< typename L, typename R, typename LB>
			struct oth_left {
				typedef typename Creator< LB, OpPower, typename Creator<Numeric, OpTimes, Numeric>::type >::type type;
				static type create(DPARM(L) l, DPARM(R) r) {
					return pow(l.l, l.r * r);
				}
			};

			template<typename L, typename R, typename LB, typename RB>
			struct apply<Expression<LB, OpPower, RB>, Numeric, L, R> : mpl::if_<
				is_same< typename remove_all< RB >::type, Numeric >,
				oth_left< L, R, LB >,
				Private::DefaultExpression< L, OpPower, R, LB, RB >
			>::type {};
		};

		template<typename WRT, typename L, typename R>
		static typename DerivType::template apply<L, R, WRT>::type deriv(WRT const& wrt, DPARM(L) left, DPARM(R) right) {
			typedef typename SKT::derivatives::Creator<L, SKT::derivatives::OpLog, double>::type log_type;
			return (log(left) * DGD(R, right, WRT, wrt) + DGD(log_type, log(left), WRT, wrt) * right) * pow(left, right);
		}

		template<typename L, typename R>
		static void print(L const& l, R const& r) {
			printf("("); l.print(); printf(" ^ "); r.print(); printf(")");
		}
	};

// /
	struct OpDivide {
		template<typename L, typename R>
		static double eval(DPARM(L) left, DPARM(R) right) {
			return left.eval() / right.eval();
		}

		struct DerivType {
			template<typename L, typename R, typename WRT>
			struct apply : Creator<
				typename Creator<
					typename Creator<R, OpTimes, typename Derivator< L, WRT >::type >::type,
					OpMinus,
					typename Creator<L, OpTimes, typename Derivator< R, WRT >::type >::type
				>::type,
				OpDivide,
				typename Creator< R, OpPower, Numeric >::type
			> {};
		};

		struct Simplify {
			SK_DERIV_SIMPLIFY_BASE(OpDivide);

			template<typename RB, typename L, typename R>
			struct apply<Zero, RB, L, R> {
				typedef Zero type;
				static type create(DPARM(L) l, DPARM(R) r) { return Zero(); }
			};

			template<typename LB, typename L, typename R>
			struct apply<LB, One, L, R> {
				typedef typename type_chooser< L >::held_type type;
				static type create(DPARM(L) l, DPARM(R) r) { return l; }
			};

			template<typename L, typename R>
			struct apply<Zero, One, L, R> {
				typedef Zero type;
				static type create(DPARM(L) l, DPARM(R) r) { return Zero(); }
			};

			template<typename LB, typename L, typename R>
			struct apply<LB, Zero, L, R> {
				BOOST_STATIC_ASSERT_MSG(false, "attempted division by Zero");
				//typedef Zero type;
				//static type create(DPARM(L) l, DPARM(R) r) { return Zero(); }
			};
		};

		template<typename WRT, typename L, typename R>
		static typename DerivType::template apply<L, R, WRT>::type deriv(WRT const& wrt, L const& left, R const& right) {
			return (right * DGD(L, left, WRT, wrt) - left * DGD(R, right, WRT, wrt)) / pow(right, Numeric(2.0));
		}

		template<typename L, typename R>
		static void print(L const& l, R const& r) {
			printf("("); l.print(); printf(" / "); r.print(); printf(")");
		}
	};

//--------------------------------------------------------------------
// proxies
//--------------------------------------------------------------------
	template<typename Target, typename WRT>
	struct DerivativeProxy<1, Target, WRT> : differentiable_class_tag {
		DerivativeProxy(Target const& target, WRT const& wrt) : target(target), wrt_(wrt) {
		}

		template<typename WRTN>
		struct GetDeravitive {
			typedef DerivativeProxy<1, Target, WRT> param_type;
			static DerivativeProxy<2, DerivativeProxy, WRTN> deriv(DPARM(param_type) dp, DPARM(WRTN) wrt) {
				return DerivativeProxy<2, DerivativeProxy, WRTN>(dp, wrt);
			}
		};

		double eval(void) const {
			return target.eval1(this->wrt1());
		};

		WRT const& wrt1(void) const {
			return wrt_;
		};

	protected:
		Target const& target;

	private:
		WRT const& wrt_;
	};


	template<typename Base, typename WRT>
	struct DerivativeProxy<2, Base, WRT> : Base {

		DerivativeProxy(Base const& base, WRT const& wrt) : Base(base), wrt_(wrt) {
		}

		template<typename WRTN>
		struct GetDeravitive {
			typedef DerivativeProxy<2, Base, WRT> param_type;
			static DerivativeProxy<3, DerivativeProxy, WRTN> deriv(DPARM(param_type) dp, DPARM(WRTN) wrt) {
				return DerivativeProxy<3, DerivativeProxy, WRTN>(dp, wrt);
			}
		};

		double eval(void) const {
			return this->target.eval2(this->wrt1(), this->wrt2());
		};

		WRT const& wrt2(void) const {
			return wrt_;
		};

	private:
		WRT const& wrt_;
	};


	template<typename Base, typename WRT>
	struct DerivativeProxy<3, Base, WRT> : Base {

		DerivativeProxy(Base const& base, WRT const& wrt) : Base(base), wrt_(wrt) {
		}

		template<typename WRTN>
		struct GetDeravitive {
			typedef DerivativeProxy<3, Base, WRT> param_type;
			static DerivativeProxy<4, DerivativeProxy, WRTN> deriv(DPARM(param_type) dp, DPARM(WRTN) wrt) {
				return DerivativeProxy<4, DerivativeProxy, WRTN>(dp, wrt);
			}
		};

		double eval(void) const {
			return this->target.eval3(this->wrt1(), this->wrt2(), this->wrt3());
		};

		WRT const& wrt3(void) const {
			return wrt_;
		};

	private:
		WRT const& wrt_;
	};


	template<typename Base, typename WRT>
	struct DerivativeProxy<4, Base, WRT> : Base {

		DerivativeProxy(Base const& base, WRT const& wrt) : Base(base), wrt_(wrt) {}

		double eval(void) const {
			return this->target.eval4(this->wrt1(), this->wrt2(), this->wrt3(), this->wrt4());
		};

		WRT const& wrt4(void) const {
			return wrt_;
		};

	private:
		WRT const& wrt_;
	};

	namespace Private {
		template<typename T, typename WRT>
		struct DerivativeType {
			typedef DerivativeProxy<1, T, WRT> type;
		};

		template<typename T>
		struct DerivativeType<T, T> {
			typedef One type;
		};
	} // end namespace Private
} } // end namespace derivatives, end namespace SKT

//--------------------------------------------------------------------
// operator definitions
//--------------------------------------------------------------------
template<typename L, typename R>
typename SKT::derivatives::Creator<L, SKT::derivatives::OpPlus, R>::type operator+(L const& l, R const& r) {
	return SKT::derivatives::Creator<L, SKT::derivatives::OpPlus, R>::create(l, r);
}

template<typename L, typename R>
typename SKT::derivatives::Creator<L, SKT::derivatives::OpMinus, R>::type operator-(L const& l, R const& r) {
	return SKT::derivatives::Creator<L, SKT::derivatives::OpMinus, R>::create(l, r);
}

template<typename L, typename R>
typename SKT::derivatives::Creator<L, SKT::derivatives::OpTimes, R>::type operator*(L const& l, R const& r) {
	return SKT::derivatives::Creator<L, SKT::derivatives::OpTimes, R>::create(l, r);
}

template<typename L, typename R>
typename SKT::derivatives::Creator<L, SKT::derivatives::OpDivide, R>::type operator/(L const& l, R const& r) {
	return SKT::derivatives::Creator<L, SKT::derivatives::OpDivide, R>::create(l, r);
}

template<typename L, typename R>
typename SKT::derivatives::Creator<L, SKT::derivatives::OpPower, R>::type pow(L const& l, R const& r) {
	return SKT::derivatives::Creator<L, SKT::derivatives::OpPower, R>::create(l, r);
};

template<typename R>
typename SKT::derivatives::Creator<SKT::derivatives::Numeric, SKT::derivatives::OpPower, R>::type exp(R const& r) {
	return SKT::derivatives::Creator<SKT::derivatives::Numeric, SKT::derivatives::OpPower, R>::create(EXP, r);
};

template<typename L>
typename SKT::derivatives::Creator<L, SKT::derivatives::OpLog, SKT::derivatives::Numeric>::type log(L const& l) {
	return SKT::derivatives::Creator<L, SKT::derivatives::OpLog, SKT::derivatives::Numeric>::create(l, 0.0);
};

//--------------------------------------------------------------------
// derivation functions
//--------------------------------------------------------------------
template<typename WRT, typename T>
typename SKT::derivatives::Derivator<T, WRT>::type D(T const & of, WRT const & wrt) {
	return DGD(T, of, WRT, wrt);
	//DGD(D_STRUCT_T, D_STRUCT_V, D_WRT_T, D_WRT_V ) D_STRUCT_T::template GetDeravitive<D_WRT_T>::deriv(D_STRUCT_V, D_WRT_V)
	//return of.deriv(wrt);
};

template<typename WRT>
double D(double of, WRT const& wrt) {
	return 0.0;
};

//--------------------------------------------------------------------
// evaluation functions
//--------------------------------------------------------------------
template<typename T>
double E(T const &t){
	return t.eval();
}

double E(double v){
	return v;
}

//--------------------------------------------------------------------
// printing functions
//--------------------------------------------------------------------
template<typename T>
void P(T const &t) {
	t.print(); printf(" = %g\n\n", E(t));
}

void P(double t) {
	printf("%g", t); printf("\n");
}

//--------------------------------------------------------------------
// macros
//--------------------------------------------------------------------
#define SK_DERIV_DEFINE_VARIABLE(name) \
struct Numeric_##name : SKT::derivatives::Numeric { \
	explicit Numeric_##name(Numeric_##name const& other) : SKT::derivatives::Numeric(other) {} \
	explicit Numeric_##name(double t) : SKT::derivatives::Numeric(t) {} \
	template<typename T> struct DerivativeType : SKT::derivatives::DerivativeType< Numeric_##name, T > {}; \
	template<typename T> struct GetDeravitive { \
		static typename DerivativeType< T >::type deriv(DPARMET(Numeric_##name), DPARMET(T)) { \
			return DerivativeType< T >::type(); \
		} \
	}; \
	void print(void) const { printf(#name"["); SKT::derivatives::Numeric::print(); printf("]"); } \
}; \
typedef Numeric_##name typeof##name

#define SK_DERIV_DECLARE_VARIABLE(name) Numeric_##name name

#define SK_DERIV_INIT_VARIABLE(name, val) name(val)

#define SK_DERIV_MEMBER_VARIABLE(name) SK_DERIV_DEFINE_VARIABLE(name); SK_DERIV_DECLARE_VARIABLE(name)

#define SK_DERIV_GLOBAL_VARIABLE(name, val) SK_DERIV_DEFINE_VARIABLE(name); Numeric_##name name(val)

#define SK_DERIV_DECLARE_DERIVATIVES(classname, expression) \
	template<typename WRT> \
	struct DerivativeType : SKT::derivatives::Private::DerivativeType<classname, WRT> {}; \
\
	template<typename WRT> struct GetDeravitive { \
		typedef typename DerivativeType< WRT >::type type; \
		static type deriv(DPARMET(classname) t, DPARMET(WRT) wrt) { return type(t, wrt); } \
	}; \
\
	double eval(void) const { \
		return E(expression); \
	} \
\
	template<typename WRT> \
	double eval1(WRT const& wrt) const { \
		return E(D(expression, wrt)); \
	} \
\
	template<typename WRT1, typename WRT2> \
	double eval2(WRT1 const& wrt1, WRT2 const& wrt2) const { \
		return E(D(D(expression, wrt1), wrt2)); \
	} \
\
	template<typename WRT1, typename WRT2, typename WRT3> \
	double eval3(WRT1 const& wrt1, WRT2 const& wrt2, WRT3 const& wrt3) const { \
		return E(D(D(D(expression, wrt1), wrt2), wrt3)); \
	} \
\
	template<typename WRT1, typename WRT2, typename WRT3, typename WRT4> \
	double eval4(WRT1 const& wrt1, WRT2 const& wrt2, WRT3 const& wrt3, WRT4 const& wrt4) const { \
		return E(D(D(D(D(expression, wrt1), wrt2), wrt3), wrt4)); \
	}

#endif //SK_DERIVIITIVES_HEADER_INCLUDED_298YWDOIKJHWOTUYHWOKJGNLIUHWT








//\
//	/*template<typename WRT>*/ \
//	/*SKT::derivatives::DerivativeProxy<1, classname, WRT> deriv(DPARMET(WRT) wrt) const { */\
//	/*	return SKT::derivatives::DerivativeProxy<1, classname, WRT>(*this, wrt); */\
//	/*}*/ \
//\
//	/*SKT::derivatives::One deriv(classname const&) const {*/ \
//	/*	return SKT::derivatives::One();*/ \
//	/*}*/