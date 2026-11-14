/**
 * @file test_gsl.cpp
 * @author Zhikui Guo (zguo@geomar.de)
 * @brief This test example comes from GSL release 2.7 chapter 38.8
 * \f{equation}
 * \begin{matrix}
 * f_1(x,y) & =  & a(1-x) \\
 * f_2(x,y) & = & b(y-x^2)
 * \end{matrix}, (a=1, b=10)
 * \f}
 * @version 0.1
 * @date 2021-12-27
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
struct rparams
{
    double a;
    double b; 
};
int rosenbrock_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  double a = ((struct rparams *) params)->a;
  double b = ((struct rparams *) params)->b;
  const double x0 = gsl_vector_get (x, 0);
  const double x1 = gsl_vector_get (x, 1);
  const double y0 = a * (1 - x0);
  const double y1 = b * (x1 - x0 * x0);
  gsl_vector_set (f, 0, y0);
  gsl_vector_set (f, 1, y1);
    return GSL_SUCCESS; 
}
int rosenbrock_df (const gsl_vector * x, void *params, gsl_matrix * J)
{
  const double a = ((struct rparams *) params)->a;
  const double b = ((struct rparams *) params)->b;

  const double x0 = gsl_vector_get (x, 0);

  const double df00 = -a;
  const double df01 = 0;
  const double df10 = -2 * b  * x0;
  const double df11 = b;

  gsl_matrix_set (J, 0, 0, df00);
  gsl_matrix_set (J, 0, 1, df01);
  gsl_matrix_set (J, 1, 0, df10);
  gsl_matrix_set (J, 1, 1, df11);

  return GSL_SUCCESS;
}
int rosenbrock_fdf (const gsl_vector * x, void *params,
                gsl_vector * f, gsl_matrix * J)
{
  rosenbrock_f (x, params, f);
  rosenbrock_df (x, params, J);

  return GSL_SUCCESS;
}
void print_state (size_t iter, gsl_multiroot_fsolver * s)
{
  printf ("iter = %3lu x = % 15.8f % 15.8f  f(x) = % .3e % .3e\n",iter,
          gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1),
          gsl_vector_get (s->f, 0), gsl_vector_get (s->f, 1));
}
void print_state (size_t iter, gsl_multiroot_fdfsolver * s)
{
  printf ("iter = %3lu x = % 15.8f % 15.8f  f(x) = % .3e % .3e\n",iter,
          gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1),
          gsl_vector_get (s->f, 0), gsl_vector_get (s->f, 1));
}
void test_gsl_release2_7_chap38_8()
{
    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;
    int status;
    size_t i, iter = 0;
    const size_t n = 2;
    struct rparams p = {1.0, 10.0};
    gsl_multiroot_function f = {&rosenbrock_f, n, &p};
    double x_init[2] = {-10, -5}; //initial guess
    gsl_vector *x = gsl_vector_alloc (n);
    gsl_vector_set (x, 0, x_init[0]);
    gsl_vector_set (x, 1, x_init[1]);
    T = gsl_multiroot_fsolver_hybrids;
    s = gsl_multiroot_fsolver_alloc (T, 2);
    gsl_multiroot_fsolver_set (s, &f, x);
    print_state (iter, s);

    do 
    {
        iter++;
        status = gsl_multiroot_fsolver_iterate (s);
        print_state (iter, s);
        if (status) /* check if solver is stuck */
        break;
        status =
        gsl_multiroot_test_residual (s->f, 1e-7);
    }while (status == GSL_CONTINUE && iter < 1000);
    printf ("status = %s\n", gsl_strerror (status));
    gsl_multiroot_fsolver_free (s); gsl_vector_free (x);
}
void test_gsl_release2_7_chap38_8_fdf()
{
    const gsl_multiroot_fdfsolver_type *T;
    gsl_multiroot_fdfsolver *s;
    int status;
    size_t i, iter = 0;
    const size_t n = 2;
    struct rparams p = {1.0, 10.0};
    gsl_multiroot_function_fdf f = {&rosenbrock_f, &rosenbrock_df, &rosenbrock_fdf, n, &p};
    double x_init[2] = {-10, -5}; //initial guess
    gsl_vector *x = gsl_vector_alloc (n);
    gsl_vector_set (x, 0, x_init[0]);
    gsl_vector_set (x, 1, x_init[1]);
    T = gsl_multiroot_fdfsolver_gnewton;
    s = gsl_multiroot_fdfsolver_alloc (T, 2);
    gsl_multiroot_fdfsolver_set (s, &f, x);
    print_state (iter, s);

    do 
    {
        iter++;
        status = gsl_multiroot_fdfsolver_iterate (s);
        print_state (iter, s);
        if (status) /* check if solver is stuck */
        break;
        status =
        gsl_multiroot_test_residual (s->f, 1e-7);
    }while (status == GSL_CONTINUE && iter < 1000);
    printf ("status = %s\n", gsl_strerror (status));
    gsl_multiroot_fdfsolver_free (s); gsl_vector_free (x);
}

/**
 * @brief See https://nlopt.readthedocs.io/en/latest/NLopt_Tutorial/
 * \f$ f=(2x + 0)^3 - (-x + 1)^3 = 0 \f$
 */
struct params
{
    double a, b, c, d;
};
int func(const gsl_vector * x, void *params, gsl_vector* f)
{
    double a = ((struct params *) params)->a;
    double b = ((struct params *) params)->b;
    double c = ((struct params *) params)->c;
    double d = ((struct params *) params)->d;
    const double x0 = gsl_vector_get (x, 0);
    const double y0 = pow(a*x0 + b, 3.0) - pow(c*x0 + d, 3.0);
    gsl_vector_set (f, 0, y0);
    return GSL_SUCCESS; 
}
int func_df (const gsl_vector * x, void *params, gsl_matrix * J)
{
  const double a = ((struct rparams *) params)->a;
  const double b = ((struct rparams *) params)->b;

  const double x0 = gsl_vector_get (x, 0);

  const double df00 = 6*(2*x0)*(2*x0) + 3*(-x0 + 1)*(-x0+1);

  gsl_matrix_set (J, 0, 0, df00);

  return GSL_SUCCESS;
}
int func_fdf (const gsl_vector * x, void *params,
                gsl_vector * f, gsl_matrix * J)
{
  func (x, params, f);
  func_df (x, params, J);

  return GSL_SUCCESS;
}
void print_state2 (size_t iter, gsl_multiroot_fsolver * s)
{
  printf ("iter = %3lu x = % 15.8f f(x) = % .3e \n",iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->f, 0));
}
void print_state2 (size_t iter, gsl_multiroot_fdfsolver * s)
{
  printf ("iter = %3lu x = % 15.8f f(x) = % .3e \n",iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->f, 0));
}
void test_NLopt_tutorial()
{
    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;
    int status;
    size_t i, iter = 0;
    const size_t n = 1;
    struct params p = {2, 0, -1, 1};
    gsl_multiroot_function f = {&func, n, &p};
    double x_init[n] = {3}; //initial guess
    gsl_vector *x = gsl_vector_alloc (n);
    gsl_vector_set (x, 0, x_init[0]);
    T = gsl_multiroot_fsolver_hybrids;
    s = gsl_multiroot_fsolver_alloc (T, n);
    gsl_multiroot_fsolver_set (s, &f, x);
    print_state2 (iter, s);

    do 
    {
        iter++;
        status = gsl_multiroot_fsolver_iterate (s);
        print_state2 (iter, s);
        if (status) /* check if solver is stuck */
        break;
        status =
        gsl_multiroot_test_residual (s->f, 1e-15);
    }while (status == GSL_CONTINUE && iter < 1000);
    printf ("status = %s\n", gsl_strerror (status));
    gsl_multiroot_fsolver_free (s); gsl_vector_free (x);
}
void test_NLopt_tutorial_fdf()
{
    const gsl_multiroot_fdfsolver_type *T;
    gsl_multiroot_fdfsolver *s;
    int status;
    size_t i, iter = 0;
    const size_t n = 1;
    struct params p = {2, 0, -1, 1};
    gsl_multiroot_function_fdf f = {&func, &func_df, &func_fdf, n, &p};
    double x_init[n] = {3}; //initial guess
    gsl_vector *x = gsl_vector_alloc (n);
    gsl_vector_set (x, 0, x_init[0]);
    T = gsl_multiroot_fdfsolver_gnewton;
    s = gsl_multiroot_fdfsolver_alloc (T, n);
    gsl_multiroot_fdfsolver_set (s, &f, x);
    print_state2 (iter, s);

    do 
    {
        iter++;
        status = gsl_multiroot_fdfsolver_iterate (s);
        print_state2 (iter, s);
        if (status) /* check if solver is stuck */
        break;
        status =
        gsl_multiroot_test_residual (s->f, 1e-15);
    }while (status == GSL_CONTINUE && iter < 1000);
    printf ("status = %s\n", gsl_strerror (status));
    gsl_multiroot_fdfsolver_free (s); gsl_vector_free (x);
}

struct quadratic_params
  {
    double a, b, c;
  };

double quadratic (double x, void *params);
double quadratic_deriv (double x, void *params);
void quadratic_fdf (double x, void *params,
                    double *y, double *dy);
int test_onedimensional_root_finding();
int main()
{
    // test_gsl_release2_7_chap38_8();
    // test_gsl_release2_7_chap38_8_fdf();
    // test_NLopt_tutorial();
    // test_NLopt_tutorial_fdf();
    test_onedimensional_root_finding();
    return 0;
}
int test_onedimensional_root_finding()
{
  int status;
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r = 0, r_expected = sqrt (5.0);
  double x_lo = 0.0, x_hi = 5.0;
  gsl_function F;
  struct quadratic_params params = {1.0, 0.0, -5.0};

  F.function = &quadratic;
  F.params = &params;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);

  printf ("using %s method\n",
          gsl_root_fsolver_name (s));

  printf ("%5s [%9s, %9s] %9s %10s %9s\n",
          "iter", "lower", "upper", "root",
          "err", "err(est)");

  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      r = gsl_root_fsolver_root (s);
      x_lo = gsl_root_fsolver_x_lower (s);
      x_hi = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lo, x_hi,
                                       0, 0.001);

      if (status == GSL_SUCCESS)
        printf ("Converged:\n");

      printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
              iter, x_lo, x_hi,
              r, r - r_expected,
              x_hi - x_lo);
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fsolver_free (s);

  return status;
}
double
quadratic (double x, void *params)
{
  struct quadratic_params *p
    = (struct quadratic_params *) params;

  double a = p->a;
  double b = p->b;
  double c = p->c;
  // printf("x = %f, %f\n",x,(a * x + b) * x + c);
  return (a * x + b) * x + c;
}

double
quadratic_deriv (double x, void *params)
{
  struct quadratic_params *p
    = (struct quadratic_params *) params;

  double a = p->a;
  double b = p->b;

  return 2.0 * a * x + b;
}

void
quadratic_fdf (double x, void *params,
               double *y, double *dy)
{
  struct quadratic_params *p
    = (struct quadratic_params *) params;

  double a = p->a;
  double b = p->b;
  double c = p->c;

  *y = (a * x + b) * x + c;
  *dy = 2.0 * a * x + b;
}