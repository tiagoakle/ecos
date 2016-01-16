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

#if CONEMODE == 0
#define MTILDE (m+2*C->nsoc)
#else
#define MTILDE (m)
#endif

/* Solve with the abs of the diagonal*/
void inline LDL_abs_dsolve(idxint n, pfloat *x, pfloat * D)
{
    idxint i = 0;
    for(i;i<n;i++)
        x[i] /= fabs(D[i]);
}

void inline preconditioner_solve(idxint nK, pfloat* Pb, pfloat *Px, kkt* KKT)
{
    /* forward - diagonal - backward solves: Px holds solution */
    LDL_lsolve2(nK, Pb, KKT->L->jc, KKT->L->ir, KKT->L->pr, Px);
    LDL_abs_dsolve(nK, Px, KKT->D);
    LDL_ltsolve(nK, Px, KKT->L->jc, KKT->L->ir, KKT->L->pr);
}


//TODO DOCUMENT input output and operation permutations
/**
 * Computes the product PKP', the input is expected in
 * the output is set in Px
 */
void inline KKT_prod( idxint n,
                      idxint m,
                      idxint p,
                      idxint nK,
                      cone* C,
                      spmat* A,
                      spmat* G,
                      idxint isinit,
                      pfloat* dx, //Work
                      pfloat* dy, //Work
                      pfloat* dz, //Work
                      idxint *Pinv,
                      pfloat *Px, //Input
                      pfloat *work2,
                      pfloat *work3,
                      pfloat *Ab,  //Assign to work4, output
                      pfloat *work5,
                      pfloat *work6)
{
    idxint i, k, l, j, kk, kItRef;
    pfloat* dPx = work2;
    pfloat* e  = work3;
    pfloat* ex = work3;
    pfloat* ey = e + n;
    pfloat* ez = e + n+p;
    pfloat* Gdx = work6;
    pfloat* truez = work5;

    /* unpermute x & copy into arrays */
    unstretch(n, p, C, Pinv, Px, dx, dy, dz);

    /* compute error term */
    k=0; j=0;

    /* 1. product on dx*/

    /* ex = A'*dy + G'*dz */
    if(A) sparseMtVm(A, dy, ex, 1, 0);
    sparseMtVm(G, dz, ex, -1, 0);
    //TODO is p the size of the rows of A or cols of A'?
    for(j=0;j<p;j++)
        ex[j] = -ex[j];
    j=0;

    /* product on dy */
    if( p > 0 ){
        /* ey =  A*dx */
        sparseMV(A, dx, ey, 1, 0);
    }

    //TODO: Remove uneccesary Gdx??
    /* --> 3. ez = G*dx - V*dz_true */
    kk = 0; j=0;
    sparseMV(G, dx, Gdx, 1, 1);
    for( i=0; i<C->lpc->p; i++ ){

        ez[kk++] =  Gdx[j++];
    }

    for( l=0; l<C->nsoc; l++ ){
        for( i=0; i<C->soc[l].p; i++ ){
            ez[kk++] =  Gdx[j++];
        }
#if CONEMODE == 0
        ez[kk] = 0;
        ez[kk+1] = 0;
        k += 2;
        kk += 2;
#endif
    }
#ifdef EXPCONE
    for(l=0; l<C->nexc; l++)
        {
            for(i=0;i<3;i++)
            {
				ez[kk++] = Gdx[j++];
            }
        }
#endif
    //Multiply with the hessian
    for( i=0; i<MTILDE; i++) { truez[i] = -Px[Pinv[n+p+i]]; }
    if( isinit == 0 ){
        scale2add(truez, ez, C);
    } else {
        vadd(MTILDE, truez, ez);
    }

    //Finally permute again and assign to output
    for( i=0; i<nK; i++) { Ab[Pinv[i]] = e[i]; }
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
//TODO rename to kkt_solve
        idxint kkt_solve_min(kkt *KKT, spmat *A, spmat *G, pfloat *Pb, pfloat *dx, pfloat *dy, pfloat *dz, idxint n,
                         idxint p, idxint m, cone *C, idxint isinit, idxint nitref)
{

#if CONEMODE == 0
#define MTILDE (m+2*C->nsoc)
#else
#define MTILDE (m)
#endif

    idxint i, k, l, j, kk, kItRef;
#if (defined STATICREG) && (STATICREG > 0)
    idxint dzoffset;
#endif
    idxint *Pinv = KKT->Pinv;
    pfloat *Px = KKT->work1;
    pfloat *dPx = KKT->work2;
    pfloat *e = KKT->work3;
    pfloat *Pe = KKT->work4;
    pfloat *truez = KKT->work5;
    pfloat *Gdx = KKT->work6;
    pfloat *ex = e;
    pfloat *ey = e + n;
    pfloat *ez = e + n + p;
    pfloat bnorm = 1.0 + norminf(Pb, n + p + MTILDE);
    pfloat nex = 0;
    pfloat ney = 0;
    pfloat nez = 0;
    pfloat nerr;
    pfloat nerr_prev = (pfloat) NAN;
    pfloat error_threshold = bnorm * LINSYSACC;
    idxint nK = KKT->PKPt->n;

    /* forward - diagonal - backward solves: Px holds solution */
    LDL_lsolve2(nK, Pb, KKT->L->jc, KKT->L->ir, KKT->L->pr, Px);
    LDL_dsolve(nK, Px, KKT->D);
    LDL_ltsolve(nK, Px, KKT->L->jc, KKT->L->ir, KKT->L->pr);

    unstretch(n, p, C, Pinv, Px, dx, dy, dz);

    return kItRef;
}
