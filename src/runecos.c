/*
 * ECOS - Embedded Conic Solver.
 * Copyright (C) 2012-14 Alexander Domahidi [domahidi@control.ee.ethz.ch],
 * Automatic Control Laboratory, ETH Zurich.
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


/* main file with example of how to run ECOS */

#include <stdio.h>

#include "ecos.h"
#include "data.h"

int main(void)
{
	/*char ver[7];*/
    idxint exitflag = ECOS_FATAL;
	pwork* mywork;
#if PROFILING > 0
	double ttotal, tsolve, tsetup;
#endif
#if PROFILING > 1
    double torder, tkktcreate, ttranspose, tfactor, tkktsolve;
#endif
	
	/* set up data */	
#ifdef EXPCONE
	mywork = ECOS_setup(n, m, p, l, ncones, q, 0, Gpr, Gjc, Gir, Apr, Ajc, Air, c, h, b);
#else
	mywork = ECOS_setup(n, m, p, l, ncones, q, Gpr, Gjc, Gir, Apr, Ajc, Air, c, h, b);
#endif
 
    if( mywork != NULL ){
	
	/* solve */	
	exitflag = ECOS_solve(mywork);
    
    /* test second solve
    exitflag = ECOS_solve(mywork); */

	/* some statistics in milliseconds */
#if PROFILING > 0
	tsolve = mywork->info->tsolve         * 1000;
	tsetup = mywork->info->tsetup         * 1000;
	ttotal = tsetup + tsolve;
#endif
#if PROFILING > 1
	torder = mywork->info->torder         * 1000;
	tkktcreate = mywork->info->tkktcreate * 1000;
	ttranspose = mywork->info->ttranspose * 1000;
	tfactor = mywork->info->tfactor       * 1000;
	tkktsolve = mywork->info->tkktsolve   * 1000;
#endif

#if PRINTLEVEL > 2
#if PROFILING > 1
	printf("ECOS timings (all times in milliseconds):\n\n");
	printf("1. Setup: %7.3f (%4.1f%%)\n", tsetup,  tsetup / ttotal*100);
	printf("2. Solve: %7.3f (%4.1f%%)\n", tsolve,  tsolve / ttotal*100);
	printf("----------------------------------\n");
	printf(" Total solve time: %7.3f ms\n\n", ttotal);
	
	printf("Detailed timings in SETUP:\n");
	printf("Create transposes: %7.3f (%4.1f%%)\n", ttranspose, ttranspose / tsetup*100);
	printf("Create KKT Matrix: %7.3f (%4.1f%%)\n", tkktcreate, tkktcreate / tsetup*100);
	printf(" Compute ordering: %7.3f (%4.1f%%)\n", torder,         torder / tsetup*100);
	printf("            Other: %7.3f (%4.1f%%)\n", tsetup-torder-tkktcreate-ttranspose,         (tsetup-torder-tkktcreate-ttranspose) / tsetup*100);
	printf("\n");

	printf("Detailed timings in SOLVE:\n");
	printf("   Factorizations: %7.3f (%4.1f%% of tsolve / %4.1f%% of ttotal)\n", tfactor,     tfactor / tsolve*100, tfactor / ttotal*100);
	printf("       KKT solves: %7.3f (%4.1f%% of tsolve / %4.1f%% of ttotal)\n", tkktsolve, tkktsolve / tsolve*100, tfactor / ttotal*100);
	printf("            Other: %7.3f (%4.1f%% of tsolve / %4.1f%% of ttotal)\n", tsolve-tkktsolve-tfactor, (tsolve-tkktsolve-tfactor) / tsolve*100, (tsolve-tkktsolve-tfactor) / ttotal*100);
	printf("\n");
#endif
#endif
	
    /* clean up memory */
	ECOS_cleanup(mywork, 0);
        
    }
    
    /* test version number
    ECOS_ver(ver);
    printf("This test has been run on ECOS version %s\n", ver);
     */
	
    /* explicitly truncate exit code */
	return (int)exitflag;
}
