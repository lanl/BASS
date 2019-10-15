/*************************** stocR.cpp **********************************
* Author:        Agner Fog
* Date created:  2006
* Last modified: 2011-08-05
* Project:       BiasedUrn
* Source URL:    www.agner.org/random
*
* Description:
* Interface of non-uniform random number generators to R-language implementation.
* This file contains source code for the class StocRBase defined in stocR.h.
*
* Copyright 2006-2011 by Agner Fog. 
* GNU General Public License http://www.gnu.org/licenses/gpl.html
*****************************************************************************/

#include "stocc.h"                     // class definition

/***********************************************************************
Fatal error exit (Replaces userintf.cpp)
***********************************************************************/

void FatalError(const char * ErrorText) {
   // This function outputs an error message and aborts the program.
   error("%s", ErrorText);             // Error exit in R.DLL
}
