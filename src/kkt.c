/*
 * ECOS - Embedded Conic Solver.
 * Copyright (C) 2012-2015 A. Domahidi [domahidi@embotech.com],
 * Automatic Control Lab, ETH Zurich & embotech GmbH, Zurich, Switzerland.
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


/* The KKT module.
 * Handles all computation related to KKT matrix:
 * - updating the matrix
 * - its factorization
 * - solving for search directions
 * - etc.
 */

#include "kkt.h"
#include "ldl.h"
#include "splamm.h"
#include "ecos.h"
#include "cone.h"

#include <math.h>


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
idxint kkt_solve(kkt* KKT, spmat* A, spmat* G, pfloat* Pb, pfloat* dx, pfloat* dy, pfloat* dz, idxint n, idxint p, idxint m, cone* C, idxint isinit, idxint nitref)
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
	idxint*  Pinv = KKT->Pinv;
	pfloat*    Px = KKT->work1;
	pfloat*   dPx = KKT->work2;
	pfloat*     e = KKT->work3;
    pfloat*    Pe = KKT->work4;
    pfloat* truez = KKT->work5;
    pfloat*   Gdx = KKT->work6;
    pfloat* ex = e;
    pfloat* ey = e + n;
    pfloat* ez = e + n+p;
    pfloat bnorm = 1.0 + norminf(Pb, n+p+MTILDE);
    pfloat nex = 0;
    pfloat ney = 0;
    pfloat nez = 0;
    pfloat nerr;
    pfloat nerr_prev = (pfloat)NAN;
    pfloat error_threshold = bnorm*LINSYSACC;
    idxint nK = KKT->PKPt->n;

	/* forward - diagonal - backward solves: Px holds solution */
	LDL_lsolve2(nK, Pb, KKT->L->jc, KKT->L->ir, KKT->L->pr, Px );
	LDL_dsolve(nK, Px, KKT->D);
	LDL_ltsolve(nK, Px, KKT->L->jc, KKT->L->ir, KKT->L->pr);

#if PRINTLEVEL > 2
    if( p > 0 ){
        PRINTTEXT("\nIR: it  ||ex||   ||ey||   ||ez|| (threshold: %4.2e)\n", error_threshold);
        PRINTTEXT("    --------------------------------------------------\n");
    } else {
        PRINTTEXT("\nIR: it  ||ex||   ||ez|| (threshold: %4.2e)\n", error_threshold);
        PRINTTEXT("    -----------------------------------------\n");
    }
#endif

	/* iterative refinement */
	for( kItRef=0; kItRef <= nitref; kItRef++ ){

        /* unpermute x & copy into arrays */
        unstretch(n, p, C, Pinv, Px, dx, dy, dz);

		/* compute error term */
        k=0; j=0;

		/* 1. error on dx*/
#if (defined STATICREG) && (STATICREG > 0)
		/* ex = bx - A'*dy - G'*dz - DELTASTAT*dx */
        for( i=0; i<n; i++ ){ ex[i] = Pb[Pinv[k++]] - DELTASTAT*dx[i]; }
#else
		/* ex = bx - A'*dy - G'*dz */
		for( i=0; i<n; i++ ){ ex[i] = Pb[Pinv[k++]]; }
#endif
        if(A) sparseMtVm(A, dy, ex, -1, 0);
        sparseMtVm(G, dz, ex, -1, 0);
        nex = norminf(ex,n);

        /* error on dy */
        if( p > 0 ){
#if (defined STATICREG) && (STATICREG > 0)
			/* ey = by - A*dx + DELTASTAT*dy */
            for( i=0; i<p; i++ ){ ey[i] = Pb[Pinv[k++]] + DELTASTAT*dy[i]; }
#else
			/* ey = by - A*dx */
			for( i=0; i<p; i++ ){ ey[i] = Pb[Pinv[k++]]; }
#endif
            sparseMV(A, dx, ey, -1, 0);
            ney = norminf(ey,p);
        }


		/* --> 3. ez = bz - G*dx + V*dz_true */
        kk = 0; j=0;
#if (defined STATICREG) && (STATICREG > 0)
		dzoffset=0;
#endif
        sparseMV(G, dx, Gdx, 1, 1);
        for( i=0; i<C->lpc->p; i++ ){
#if (defined STATICREG) && (STATICREG > 0)
            ez[kk++] = Pb[Pinv[k++]] - Gdx[j++] + DELTASTAT*dz[dzoffset++];
#else
			ez[kk++] = Pb[Pinv[k++]] - Gdx[j++];
#endif
        }
        for( l=0; l<C->nsoc; l++ ){
            for( i=0; i<C->soc[l].p; i++ ){
#if (defined STATICREG) && (STATICREG > 0)
                ez[kk++] = i<(C->soc[l].p-1) ? Pb[Pinv[k++]] - Gdx[j++] + DELTASTAT*dz[dzoffset++] : Pb[Pinv[k++]] - Gdx[j++] - DELTASTAT*dz[dzoffset++];
#else
				ez[kk++] = Pb[Pinv[k++]] - Gdx[j++];
#endif
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
#if (defined STATICREG) && (STATICREG > 0)
                ez[kk++] = Pb[Pinv[k++]] - Gdx[j++] + DELTASTAT*dz[dzoffset++];
#else
				ez[kk++] = Pb[Pinv[k++]] - Gdx[j++];
#endif
            }
        }
#endif
        for( i=0; i<MTILDE; i++) { truez[i] = Px[Pinv[n+p+i]]; }
        if( isinit == 0 ){
            scale2add(truez, ez, C);
        } else {
            vadd(MTILDE, truez, ez);
        }
        nez = norminf(ez,MTILDE);


#if PRINTLEVEL > 2
        if( p > 0 ){
            PRINTTEXT("    %2d  %3.1e  %3.1e  %3.1e\n", (int)kItRef, nex, ney, nez);
        } else {
            PRINTTEXT("    %2d  %3.1e  %3.1e\n", (int)kItRef, nex, nez);
        }
#endif

        /* maximum error (infinity norm of e) */
        nerr = MAX( nex, nez);
        if( p > 0 ){ nerr = MAX( nerr, ney ); }

        /* CHECK WHETHER REFINEMENT BROUGHT DECREASE - if not undo and quit! */
        if( kItRef > 0 && nerr > nerr_prev ){
            /* undo refinement */
            for( i=0; i<nK; i++ ){ Px[i] -= dPx[i]; }
            kItRef--;
            break;
        }

        /* CHECK WHETHER TO REFINE AGAIN */
        if( kItRef == nitref || ( nerr < error_threshold ) || ( kItRef > 0 && nerr_prev < IRERRFACT*nerr ) ){
            break;
        }
        nerr_prev = nerr;

        /* permute */
        for( i=0; i<nK; i++) { Pe[Pinv[i]] = e[i]; }

        /* forward - diagonal - backward solves: dPx holds solution */
        LDL_lsolve2(nK, Pe, KKT->L->jc, KKT->L->ir, KKT->L->pr, dPx);
        LDL_dsolve(nK, dPx, KKT->D);
        LDL_ltsolve(nK, dPx, KKT->L->jc, KKT->L->ir, KKT->L->pr);

        /* add refinement to Px */
        for( i=0; i<nK; i++ ){ Px[i] += dPx[i]; }
	}

#if PRINTLEVEL > 2
    PRINTTEXT("\n");
#endif

	/* copy solution out into the different arrays, permutation included */
	unstretch(n, p, C, Pinv, Px, dx, dy, dz);

    return kItRef;
}