#ifndef GADGETRON_CUSTOM_EXPORT_H_
#define GADGETRON_CUSTOM_EXPORT_H_

#if defined (WIN32)
#if defined (__BUILD_GADGETRON_CUSTOM__)
#define EXPORTGADGETS_CUSTOM __declspec(dllexport)
#else
#define EXPORTGADGETS_CUSTOM __declspec(dllimport)
#endif
#else
#define EXPORTGADGETS_CUSTOM
#endif

#endif /* GADGETRON_CUSTOM_EXPORT_H_ */
