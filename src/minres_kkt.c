/*
 * 2016 Santiago Akle [tiagoakle@gmail.com],
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "kkt.h"
#include "ldl.h"
#include "splamm.h"
#include "ecos.h"
#include "cone.h"

#include <math.h>
#include <float.h>

#if CONEMODE == 0
#define MTILDE (m+2*C->nsoc)
#else
#define MTILDE (m)
#endif

//TODO: Modify these depending on the compiler precision
# define EPSILON DBL_EPSILON
# define REAL_MAX DBL_MAX


//AHEAD DECLARATION
//TODO move to minres.h
void minres(idxint nK, pfloat *x, pfloat *b, kkt *KKT, spmat *A, spmat *G, cone *C, idxint isinit,
            idxint n_out,
            idxint m,
            idxint p,
            idxint itnlimit,
            pfloat rtol,
            idxint show);

/* Solve with the abs of the diagonal*/
void LDL_abs_dsolve(idxint n, pfloat *x, pfloat *D) {
    idxint i = 0;
    for (i = 0; i < n; i++)
        x[i] /= fabs(D[i]);
}

void preconditioner_solve(idxint nK, pfloat *Pb, pfloat *Px, kkt *KKT) {
    /* forward - diagonal - backward solves: Px holds solution */
    LDL_lsolve2(nK, Pb, KKT->L->jc, KKT->L->ir, KKT->L->pr, Px);
    LDL_abs_dsolve(nK, Px, KKT->D);
    LDL_ltsolve(nK, Px, KKT->L->jc, KKT->L->ir, KKT->L->pr);
}


//TODO DOCUMENT input output and operation permutations
/**
 * Computes the product PKP', where P is the amd permutation,
 * and K is the sparsified KKT matrix.
 * The input is expected in Px
 * the output is set in y.
 *
 * The following work vectors are necesary
 *
 */

//Requires a work vector of size n+m+p
//A result vector of size n+m+p+2*n_soc_cones
void KKT_prod(idxint n,
              idxint m,
              idxint p,
              idxint nK,
              cone *C,
              spmat *A,
              spmat *G,
              pfloat *Px,    //Input size nK
              pfloat *y,     //Output size nK
              idxint isinit, //Multiply by hessian or identity 1 for hessian 0 for identity
              idxint *Pinv,  //Permutation vector
              pfloat *work1, //work vector of size nK
              pfloat *work2,  //work vector of size nK,
              pfloat *work3,  //work vector of size m
              pfloat *work4)  //work vector of size m+2*C->nsoc (MTILDE)

{
    idxint i, k, l, j, kk, kItRef;

    //Variables to store the result of the permutation
    pfloat *dx = work1;
    pfloat *dy = work1 + n;
    pfloat *dz = work1 + n + p;

    //Variables to store the result of the product with the matrices
    pfloat *e = work2;
    pfloat *ex = e;
    pfloat *ey = e + n;
    pfloat *ez = e + n + p;

    pfloat *Gdx = work3;
    pfloat *truez = work4;

    /* unpermute x & copy into arrays */
    unstretch(n, p, C, Pinv, Px, dx, dy, dz);

    /* compute error term */
    k = 0;
    j = 0;

    /* 1. product on dx*/

    /* ex = A'*dy + G'*dz */
    if (A) sparseMtVm(A, dy, ex, 1, 0);
    sparseMtVm(G, dz, ex, -1, 0);
    //We need to change the sign because spareseMtVm produces -A'dy
    for (j = 0; j < n; j++)
        ex[j] = -ex[j];

    /* product on dy */
    if (p > 0) {
        /* ey =  A*dx */
        sparseMV(A, dx, ey, 1, 1);
    }

    //TODO: Remove uneccesary Gdx??
    /* --> 3. ez = G*dx - V*dz_true */
    kk = 0;
    j = 0;
    sparseMV(G, dx, Gdx, 1, 1);
    for (i = 0; i < C->lpc->p; i++) {

        ez[kk++] = Gdx[j++];
    }

    for (l = 0; l < C->nsoc; l++) {
        for (i = 0; i < C->soc[l].p; i++) {
            ez[kk++] = Gdx[j++];
        }
#if CONEMODE == 0
        ez[kk] = 0;
        ez[kk + 1] = 0;
        k += 2;
        kk += 2;
#endif
    }
#ifdef EXPCONE
    for (l = 0; l < C->nexc; l++) {
        for (i = 0; i < 3; i++) {
            ez[kk++] = Gdx[j++];
        }
    }
#endif
    //Multiply with the hessian
    for (i = 0; i < MTILDE; i++) { truez[i] = -Px[Pinv[n + p + i]]; }
    if (isinit == 0) {
        scale2add(truez, ez, C);
    } else {
        vadd(MTILDE, truez, ez);
    }

    //Finally permute again and assign to output
    for (i = 0; i < nK; i++) { y[Pinv[i]] = e[i]; }
}

/**
 * Solves the permuted KKT system and returns the unpermuted search directions.
 *
 * On entry, the factorization of the permuted KKT matrix, PKPt,
 * is assumed to be up to date (call kkt_factor beforehand to achieve this).
 * The right hand side, Pb, is assumed to be already permuted.
 *
 * On exit, the resulting search directions are written into dx, dy and dz,
 * where these variables are permuted back to the original ordering.
 *
 * KKT->nitref iterative refinement steps are applied to solve the linear system.
 *
 * Returns the number of iterative refinement steps really taken.
 */
idxint kkt_solve(kkt *KKT, spmat *A, spmat *G, pfloat *Pb, pfloat *dx, pfloat *dy, pfloat *dz, idxint n,
                 idxint p, idxint m, cone *C, idxint isinit, idxint nitref) {

#if CONEMODE == 0
#define MTILDE (m+2*C->nsoc)
#else
#define MTILDE (m)
#endif

    idxint i, k, l, j, kk, kItRef;

    idxint *Pinv = KKT->Pinv;
    pfloat *Px = KKT->work13; //Minres output
    idxint nK = KKT->PKPt->n;

    /* forward - diagonal - backward solves: Px holds solution */
    //LDL_lsolve2(nK, Pb, KKT->L->jc, KKT->L->ir, KKT->L->pr, Px);
    //LDL_dsolve(nK, Px, KKT->D);
    //LDL_ltsolve(nK, Px, KKT->L->jc, KKT->L->ir, KKT->L->pr);

    idxint show = 1;
    minres(nK, Px, Pb, KKT, A, G, C, isinit, n, m, p, 20, 1E-15, show);

    //Product against K to verify the quality of the solution
    pfloat* prod_output = KKT->work5;
    pfloat* minres_sol = KKT->work6;
    //Move the original input to work6
    for(i=0;i<nK;i++)
        minres_sol[i] = Px[i];

    KKT_prod(n, m, p, nK, C, A, G, minres_sol, prod_output, isinit, Pinv, KKT->work1, KKT->work2, KKT->work3, KKT->work4);

    //for(i=0 ;i< nK ;i++) PRINTTEXT("Ax_minres[i]: %g b[i] %g\n",prod_output[i],Pb[i]);

    vsubscale(nK, 1.0, Pb, prod_output);
    pfloat norm_res = norm2(prod_output, nK);
    PRINTTEXT("Norm residual after minres %g\n", norm_res);

    unstretch(n, p, C, Pinv, minres_sol, dx, dy, dz);
    kItRef = 0;

    return kItRef;
}

void minres(idxint nK,
            pfloat *x,
            pfloat *b,
            kkt *KKT,
            spmat *A,
            spmat *G,
            cone *C,
            idxint isinit,
            idxint n_out,
            idxint m,
            idxint p,
            idxint itnlimit,
            pfloat rtol,
            idxint show) {

    //Assign the working memory
    pfloat *y = KKT->work6;
    pfloat *r1 = KKT->work7;
    pfloat *r2 = KKT->work8;
    pfloat *w = KKT->work9;
    pfloat *w1 = KKT->work10;
    pfloat *w2 = KKT->work11;
    pfloat *v = KKT->work12;

    idxint i;
    idxint n = nK;
    idxint istop, itn, done;
    istop = 0;
    itn = 0;
    done = 0;
    pfloat rnorm  = 0;
    pfloat ynorm  = 0;
    pfloat Anorm  = 0;
    pfloat Arnorm = 0;
    pfloat Acond  = 0;
    pfloat beta1  = 0;

    //Initialize the output to zeros
    for (i = 0; i < n; i++) {
        x[i] = 0;
    }

    // Set up y and v for the first Lanczos vector v1
    for (i = 0; i < n; i++) {
        r1[i] = b[i];
    }

    //Solve the first preconditioned iteration
    preconditioner_solve(n, b, y, KKT);
    //beta1 = b'y
    beta1 = eddot(n, b, y);

    if (beta1 <= 0) {
        if (beta1 < 0) istop = 9;
        show = 1;
        done = 1;
    }

    //Normalize y to get v1 later
    if (beta1 > 0) beta1 = sqrt(beta1);

    //Initialize other quantitites

    pfloat oldeps, delta, gbar, root, gamma, denom, epsa, epsx, epsr, diag, test1, test2, prnt;
    pfloat t1, t2;

    pfloat s = 0;
    pfloat z = 0;
    pfloat alfa = 0;

    pfloat oldb = 0;
    pfloat beta = beta1;
    pfloat dbar = 0;
    pfloat epsln = 0;
    pfloat qrnorm = beta1;
    pfloat phibar = beta1;
    pfloat rhs1 = beta1;
    pfloat rhs2 = 0;
    pfloat tnorm2 = 0;
    pfloat gmax = 0;
    pfloat gmin = REAL_MAX;
    pfloat cs = -1;
    pfloat sn = 0;
    pfloat phi;

    for (i = 0; i < n; i++) {
        w[i]  = 0;
        w2[i] = 0;
        r2[i] = r1[i];
    }


    //Main iteration loop
    if (done == 0) {
        while (itn < itnlimit) {
            itn += 1;
            //scale v<-1/beta y
            s = 1 / beta;
            for (i = 0; i < n; i++)
                v[i] = s * y[i];

            //KKT prod y<- Av
            KKT_prod(n_out, m, p, nK, C, A, G, v, y, isinit, KKT->Pinv, KKT->work1, KKT->work2, KKT->work3, KKT->work4);
            if (itn >= 2) {
                vsubscale(n, beta / oldb, r1, y);
            }

            alfa = eddot(n, v, y);
            vsubscale(n, alfa / beta, r2, y);
            //r1<-r2<-y
            for (i = 0; i < n; i++) {
                r1[i] = r2[i];
                r2[i] = y[i];
            }
            //Preconditioned solve
            preconditioner_solve(n, r2, y, KKT);
            oldb = beta;
            beta = eddot(n, r2, y);
            if (beta < 0) {
                istop = 9;
                break;
            }
            beta = sqrt(beta);
            tnorm2 = tnorm2 + alfa * alfa + oldb * oldb + beta * beta;

            if (itn == 1) {
                if (beta / beta1 <= 10 * EPS) //Then beta2 == 0 or close to it
                    istop = -1;
            }

            //Apply the previous rotation Qk-1
            oldeps = epsln;
            delta = cs * dbar + sn*alfa;
            gbar = sn * dbar - cs * alfa;
            epsln = sn * beta;
            dbar = -cs * beta;
            root = sqrt(gbar * gbar + dbar * dbar);
            Arnorm = phibar * root;

            //Compute the plane rotation Qk

            gamma = sqrt(gbar * gbar + beta * beta);
            gamma = fmax(gamma, EPSILON);
            cs = gbar / gamma;
            sn = beta / gamma;
            phi = cs * phibar;
            phibar = sn * phibar;

            //Update x
            denom = 1 / gamma;
            //w1<-w2<-w
            for (i = 0; i < n; i++) {
                w1[i] = w2[i];
                w2[i] = w[i];
                w[i] = (v[i] - oldeps * w1[i] - delta * w2[i]) * denom;
                x[i] = x[i] + phi * w[i];
            }

            //Go round again
            gmax = fmax(gmax, gamma);
            gmin = fmin(gmin, gamma);
            z = rhs1 / gamma;
            rhs1 = rhs2 - delta * z;
            rhs2 = -epsln * z;

            //Estimate various norms
            Anorm = sqrt(tnorm2);
            ynorm = norm2(x, n);
            epsa = Anorm * EPSILON;
            epsx = Anorm * ynorm * EPSILON;
            epsr = Anorm * ynorm * rtol;
            diag = gbar;

            if (diag == 0) diag = epsa;

            qrnorm = phibar;
            rnorm = qrnorm;
            test1 = rnorm / (Anorm * ynorm);
            test2 = root / Anorm;

            //Estimate cond A
            Acond = gmax / gmin;

            //Verify stopping criteria
            if (istop == 0) //this fails rarely if KKT = const*I (and never in ecos)
            {
                t1 = 1 + test1;
                t2 = 1 + test2;
                if (t2 <= 1) istop = 2;
                if (t1 <= 1) istop = 1;
                if (itn >= itnlimit) istop = 6;
                if (Acond >= 0.1 / EPSILON) istop = 4;
                if (epsx >= beta1) istop = 3;
                if (test2 <= rtol) istop = 2;
                if (test1 <= rtol) istop = 1;
            }

            prnt = 0;
            if (n <= 40) prnt = 1;
            if (itn <= 10) prnt = 1;
            if (itn >= itnlimit - 10) prnt = 1;
            if (itn % 10 == 0) prnt = 1;
            if (qrnorm <= 10 * epsx) prnt = 1;
            if (qrnorm <= 10 * epsr) prnt = 1;
            if (Acond <= 1E-2 / EPSILON) prnt = 1;
            if (istop != 0) prnt = 1;

            //Print and debug code
            if (prnt != 0 & show != 0) {
                if (itn % 10 == 0) PRINTTEXT(" ");
                PRINTTEXT("%6ld %12.5e %10.3e %10.3e %8.1e %8.1e %8.1e\n", itn,
                          x[0],
                          test1,
                          test2,
                          Anorm,
                          Acond,
                          gbar / Anorm);
            }
            //End of print and debug code

            if (istop != 0) break;
        }//End main loop
    }//end if done==0 from above while
    //Print final status
    if(show)
        PRINTTEXT("Final Iteration %ld: istop = %ld\n", itn, istop);
}