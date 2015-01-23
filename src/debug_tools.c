#include "debug_tools.h"

void allocateHistory(pwork *w)
{
	w->info->hist_presy = (pfloat*)MALLOC(w->stgs->maxit(sizeof(pfloat)));	
	w->info->hist_presz = (pfloat*)MALLOC(w->stgs->maxit(sizeof(pfloat)));
    w->info->hist_dres = (pfloat*)MALLOC(w->stgs->maxit(sizeof(pfloat)));
    w->info->hist_gres= (pfloat*)MALLOC(w->stgs->maxit(sizeof(pfloat)));
    w->info->hist_mu= (pfloat*)MALLOC(w->stgs->maxit(sizeof(pfloat)));
    w->info->hist_expmu= (pfloat*)MALLOC(w->stgs->maxit(sizeof(pfloat)));
    w->info->hist_sigma= (pfloat*)MALLOC(w->stgs->maxit(sizeof(pfloat)));
    w->info->hist_tau= (pfloat*)MALLOC(w->stgs->maxit(sizeof(pfloat)));
    w->info->hist_kappa= (pfloat*)MALLOC(w->stgs->maxit(sizeof(pfloat)));
    w->info->hist_hpresy= (pfloat*)MALLOC(w->stgs->maxit(sizeof(pfloat)));    
    w->info->hist_hpresz= (pfloat*)MALLOC(w->stgs->maxit(sizeof(pfloat)));
    w->info->hist_hdres= (pfloat*)MALLOC(w->stgs->maxit(sizeof(pfloat)));
}

void freeHistory(pwork *w)
{
	FREE(w->info->hist_presy);	
	FREE(w->info->hist_presz);
    FREE(w->info->hist_dres);
    FREE(w->info->hist_gres);
    FREE(w->info->hist_mu);
    FREE(w->info->hist_expmu);
    FREE(w->info->hist_sigma);
    FREE(w->info->hist_tau);
    FREE(w->info->hist_kappa);
    FREE(w->info->hist_hpresy);    
    FREE(w->info->hist_hpresz);
    FREE(w->info->hist_hdres);
}
