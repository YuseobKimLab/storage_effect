#ifndef __OPTION_H_
#define __OPTION_H_

//----------------------------------
#define _block_ 30
#define _sampling_generation_ 5000
#define _nSample_ 100
#define _PRINTOUT_EVERY_ 10000

//----------------------------------
#if defined _SBE_
    #define _MUT1_
	#define _SELSB_

#elif defined _SBC_
	#define _MUT12_
    #define _SELC_

#elif defined _SBN_
	#define _MUT12_
    #define _SELN_

#elif defined _RSE_
	#define _MUT12_
	#define _SELRS_

#elif defined _RSC_
	#define _MUT12_
	#define _SELC_

#elif defined _RSN_
	#define _MUT12_
	#define _SELN_

#endif
//----------------------------------

#define PRINT_MUTATION 0
#define PRINT_WHERE 0
#define PRINT_TML 0
#define PRINT_NONSYNSITE 0
#define PRINT_DELTA 0
#define PRINT_SELOPT 0
#define PRINT_DIST_SELOPT 0
#define PRINT_RECOM 0
#define PRINT_SELECTION 0
//----------------------------------


#endif
