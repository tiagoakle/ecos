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

#include "kkt.h"
#include "ldl.h"
#include "splamm.h"
#include "ecos.h"
#include "cone.h"

/* Factorization of KKT matrix. Just a wrapper for some LDL code */
#if PROFILING > 1
idxint kkt_factor(kkt* KKT, pfloat eps, pfloat delta, pfloat *t1, pfloat* t2)
#else
idxint kkt_factor(kkt* KKT, pfloat eps, pfloat delta)
#endif
{
    idxint nd;

    /* returns n if successful, k if D (k,k) is zero */
    nd = LDL_numeric2(
            KKT->PKPt->n,	/* K and L are n-by-n, where n >= 0 */
            KKT->PKPt->jc,	/* input of size n+1, not modified */
            KKT->PKPt->ir,	/* input of size nz=Kjc[n], not modified */
            KKT->PKPt->pr,	/* input of size nz=Kjc[n], not modified */
            KKT->L->jc,		/* input of size n+1, not modified */
            KKT->Parent,	/* input of size n, not modified */
            KKT->Sign,      /* input, permuted sign vector for regularization */
            eps,            /* input, inverse permutation vector */
            delta,          /* size of dynamic regularization */
            KKT->Lnz,		/* output of size n, not defn. on input */
            KKT->L->ir,		/* output of size lnz=Lp[n], not defined on input */
            KKT->L->pr,		/* output of size lnz=Lp[n], not defined on input */
            KKT->D,			/* output of size n, not defined on input */
            KKT->work1,		/* workspace of size n, not defn. on input or output */
            KKT->Pattern,   /* workspace of size n, not defn. on input or output */
            KKT->Flag	    /* workspace of size n, not defn. on input or output */
#if PROFILING > 1
            , t1, t2
#endif
    );

    return nd == KKT->PKPt->n ? KKT_OK : KKT_PROBLEM;
}


/**
 * Updates the permuted KKT matrix by copying in the new scalings.
 */
void kkt_update(spmat* PKP, idxint* P, cone *C)
{
    idxint i, j, k, conesize;
    pfloat eta_square, *q;
#if CONEMODE == 0
    pfloat d1, u0, u1, v1;
    idxint conesize_m1;
#else
    pfloat a, w, c, d, eta_square_d, qj;
    idxint thiscolstart;
#endif

    /* LP cone */
    for( i=0; i < C->lpc->p; i++ ){ PKP->pr[P[C->lpc->kkt_idx[i]]] = -C->lpc->v[i] - DELTASTAT; }

    /* Second-order cone */
    for( i=0; i<C->nsoc; i++ ){

#if CONEMODE == 0
        getSOCDetails(&C->soc[i], &conesize, &eta_square, &d1, &u0, &u1, &v1, &q);
        conesize_m1 = conesize - 1;

        /* D */
        PKP->pr[P[C->soc[i].Didx[0]]] = -eta_square * d1 - DELTASTAT;
        for (k=1; k < conesize; k++) {
            PKP->pr[P[C->soc[i].Didx[k]]] = -eta_square - DELTASTAT;
        }

        /* v */
        j=1;
        for (k=0; k < conesize_m1; k++) {
            PKP->pr[P[C->soc[i].Didx[conesize_m1] + j++]] = -eta_square * v1 * q[k];
        }
        PKP->pr[P[C->soc[i].Didx[conesize_m1] + j++]] = -eta_square;

        /* u */
        PKP->pr[P[C->soc[i].Didx[conesize_m1] + j++]] = -eta_square * u0;
        for (k=0; k < conesize_m1; k++) {
            PKP->pr[P[C->soc[i].Didx[conesize_m1] + j++]] = -eta_square * u1 * q[k];
        }
        PKP->pr[P[C->soc[i].Didx[conesize_m1] + j++]] = +eta_square + DELTASTAT;
#endif

#if CONEMODE > 0
        conesize = C->soc[i].p;
        eta_square = C->soc[i].eta_square;
        a = C->soc[i].a;
        w = C->soc[i].w;
        c = C->soc[i].c;
        d = C->soc[i].d;
        q = C->soc[i].q;
        eta_square_d = eta_square * d;

        /* first column - only diagonal element */
        PKP->pr[P[C->soc[i].colstart[0]]] = -eta_square * (a*a + w);

        /* next conesize-1 columns */
        for (j=1; j<conesize; j++) {

            thiscolstart = C->soc[i].colstart[j];

            /* first element in column (=c*q) */
            qj = q[j-1];
            PKP->pr[P[thiscolstart]] = -eta_square * c * qj;

            /* the rest of the column (=I + d*qq') */
            for (k=1; k<j; k++) {
                PKP->pr[P[thiscolstart+k]] = -eta_square_d * q[k-1]*qj;      /* super-diagonal elements */
            }
            PKP->pr[P[thiscolstart+j]] = -eta_square * (1.0 +  d * qj*qj);   /* diagonal element */
        }
#endif
    }
#if defined EXPCONE
    /* Exponential cones */
    for( i=0; i < C->nexc; i++){
        PKP->pr[P[C->expc[i].colstart[0]]]   = -C->expc[i].v[0]-DELTASTAT;
        PKP->pr[P[C->expc[i].colstart[1]]]   = -C->expc[i].v[1];
        PKP->pr[P[C->expc[i].colstart[1]+1]] = -C->expc[i].v[2]-DELTASTAT;
        PKP->pr[P[C->expc[i].colstart[2]]]   = -C->expc[i].v[3];
        PKP->pr[P[C->expc[i].colstart[2]+1]] = -C->expc[i].v[4];
        PKP->pr[P[C->expc[i].colstart[2]+2]] = -C->expc[i].v[5]-DELTASTAT;
    }
#endif


}

/**
 * Initializes the (3,3) block of the KKT matrix to produce the matrix
 *
 * 		[0  A'  G']
 * K =  [A  0   0 ]
 *      [G  0  -I ]
 *
 * It is assumed that the A,G have been already copied in appropriately,
 * and that enough memory has been allocated (this is done in preproc.c module).
 *
 * Note that the function works on the permuted KKT matrix.
 */
void kkt_init(spmat* PKP, idxint* P, cone *C)
{
    idxint i, j, k, conesize;
    pfloat eta_square, *q;
#if CONEMODE == 0
    pfloat d1, u0, u1, v1;
    idxint conesize_m1;
#else
    pfloat a, w, c, d, eta_square_d, qj;
    idxint thiscolstart;
#endif

    /* LP cone */
    for( i=0; i < C->lpc->p; i++ ){ PKP->pr[P[C->lpc->kkt_idx[i]]] = -1.0; }

    /* Second-order cone */
    for( i=0; i<C->nsoc; i++ ){

#if CONEMODE == 0
        getSOCDetails(&C->soc[i], &conesize, &eta_square, &d1, &u0, &u1, &v1, &q);
        conesize_m1 = conesize - 1;

        /* D */
        PKP->pr[P[C->soc[i].Didx[0]]] = -1.0;
        for (k=1; k < conesize; k++) {
            PKP->pr[P[C->soc[i].Didx[k]]] = -1.0;
        }

        /* v */
        j=1;
        for (k=0; k < conesize_m1; k++) {
            PKP->pr[P[C->soc[i].Didx[conesize_m1] + j++]] = 0.0;
        }
        PKP->pr[P[C->soc[i].Didx[conesize_m1] + j++]] = -1.0;

        /* u */
        PKP->pr[P[C->soc[i].Didx[conesize_m1] + j++]] = 0.0;
        for (k=0; k < conesize_m1; k++) {
            PKP->pr[P[C->soc[i].Didx[conesize_m1] + j++]] = 0.0;
        }
        PKP->pr[P[C->soc[i].Didx[conesize_m1] + j++]] = +1.0;
#endif

#if CONEMODE > 0
        conesize = C->soc[i].p;
        eta_square = C->soc[i].eta_square;
        a = C->soc[i].a;
        w = C->soc[i].w;
        c = C->soc[i].c;
        d = C->soc[i].d;
        q = C->soc[i].q;
        eta_square_d = eta_square * d;

        /* first column - only diagonal element */
        PKP->pr[P[C->soc[i].colstart[0]]] = -1.0;

        /* next conesize-1 columns */
        for (j=1; j<conesize; j++) {

            thiscolstart = C->soc[i].colstart[j];

            /* first element in column (=c*q) */
            qj = q[j-1];
            PKP->pr[P[thiscolstart]] = 0.0;

            /* the rest of the column (=I + d*qq') */
            for (k=1; k<j; k++) {
                PKP->pr[P[thiscolstart+k]] = 0.0;      /* super-diagonal elements */
            }
            PKP->pr[P[thiscolstart+j]] = -1.0;   /* diagonal element */
        }
#endif
    }
}

