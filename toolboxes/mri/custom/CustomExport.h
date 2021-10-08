/** \file       CustomExport.h
    \brief      Implement windows export/import for custom toolbox
    \author     Souheil Inati
*/

#pragma once

#if defined (WIN32)
    #if defined (__BUILD_GADGETRON_CUSTOM__)
        #define EXPORTCUSTOM __declspec(dllexport)
    #else
        #define EXPORTCUSTOM __declspec(dllimport)
    #endif
#else
    #define EXPORTCUSTOM
#endif
