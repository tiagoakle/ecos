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

//Ahead declaration
void LDL_abs_d_solve(idxint n, pfloat *x, pfloat * D)
{
    idxint i = 0;
    for(i;i<n;i++)
        x[i] /= fabs(D[i]);
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
