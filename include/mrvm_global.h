#ifndef MRVM_GLOBAL_H
#define MRVM_GLOBAL_H

#include <QtCore/qglobal.h>

#ifdef MRVM_LIB
# define MRVM_EXPORT Q_DECL_EXPORT
#else
# define MRVM_EXPORT Q_DECL_IMPORT
#endif

#endif // MRVM_GLOBAL_H
