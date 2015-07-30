#ifndef MLTOOLS_GLOBAL_H
#define MLTOOLS_GLOBAL_H

#include <QtCore/qglobal.h>

#if defined(MLTOOLS_LIBRARY)
#  define MLTOOLSSHARED_EXPORT Q_DECL_EXPORT
#else
#  define MLTOOLSSHARED_EXPORT Q_DECL_IMPORT
#endif

#endif // MLTOOLS_GLOBAL_H
