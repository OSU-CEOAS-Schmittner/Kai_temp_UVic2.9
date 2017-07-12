solve.h
#if defined O_embm_explicit || defined O_embm_explicit_q
# if defined O_embm_adv_q
# endif
#endif
#if !defined O_embm_explicit || defined O_embm_explicit_q
# if defined O_embm_solve2x
# else
# endif
# if defined O_embm_solve2y
# else
# endif
# if defined O_embm_mgrid
# elif defined O_embm_slap
# elif defined O_embm_adi
# elif defined O_embm_essl
# elif defined O_embm_sparskit
# endif
#endif
#if defined O_embm_solve2x
#endif
#if defined O_embm_solve2y
#endif
#if defined O_embm_solve2x || defined O_embm_solve2y
#endif
