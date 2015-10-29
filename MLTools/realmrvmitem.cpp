#include "stdafx.h"
#include "mrvm.h"

RealMRVMItem::RealMRVMItem(const QString& name)
    :MRVMItem(name)
{

}

RealMRVMItem::~RealMRVMItem()
{

}

double* RealMRVMItem::values(int index)
{
    double* value = new double[1];

    value[0] = m_values[index].toDouble();
    return value;
}

void RealMRVMItem::setValues(int index, double*& value)
{
    m_values[index] = QString::number(value[0]);
    delete[] value;
    value = NULL;
}

int RealMRVMItem::columnCount()
{
    return 1;
}

MRVMItem::MRVMValueType RealMRVMItem::mRVMValueType() const
{
    return MRVMItem::Real;
}

QString RealMRVMItem::type() const
{
    return "RealMRVMItem";
}
