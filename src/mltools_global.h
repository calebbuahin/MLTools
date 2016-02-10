#ifndef MLTOOLS_GLOBAL_H
#define MLTOOLS_GLOBAL_H

#include <QtCore/qglobal.h>

#ifdef MLTOOLS_LIB
# define MLTOOLS_EXPORT Q_DECL_EXPORT
#else
# define MLTOOLS_EXPORT Q_DECL_IMPORT
#endif

#endif // MLTOOLS_GLOBAL_H
