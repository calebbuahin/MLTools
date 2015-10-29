#include "stdafx.h"
#include "mrvm.h"

RealArrayMRVMItem::RealArrayMRVMItem(const QString& name)
    :MRVMItem(name), m_columnCount(-1)
{

}

RealArrayMRVMItem::~RealArrayMRVMItem()
{

}

double* RealArrayMRVMItem::values(int index)
{
    if(index < m_values.length())
    {
        QString s = m_values[index];

        s = s.replace("[","").replace("]","");
        QStringList args = s.split(":");

        assert(args.length() == columnCount());

        if(args.length() > 0)
        {
           int length = args.length();

           double* values = new double[args.length()];

           for(int i = 0 ; i < length ; i++)
           {
              values[i] = args[i].toDouble();
           }

           return values;
        }
    }

    return NULL;
}

void RealArrayMRVMItem::setValues(int index, double*& values)
{
    if(columnCount() > 0 && index < m_values.count())
    {
        QString s = QString("%").arg(values[0]);

        if(m_columnCount > 1)
        {
            for(int i = 1; i < m_columnCount; i++)
            {
                s = s + QString(",%").arg(values[i]);
            }
        }

        s = "[" + s + "]";
        m_values[index] = s;
    }

    delete[] values;
    values = NULL;
}

int RealArrayMRVMItem::columnCount()
{
    if(m_columnCount <= -1)
    {
        if(m_values.length() > 0)
        {
            QString s = m_values[0];
            s.replace("[","").replace("]","");
            QStringList args = s.split(":");

            m_columnCount = args.length();
        }
    }

    return m_columnCount;
}

MRVMItem::MRVMValueType RealArrayMRVMItem::mRVMValueType() const
{
    return MRVMItem::Real;
}

QString RealArrayMRVMItem::type() const
{
    return "RealArrayMRVMItem";
}
