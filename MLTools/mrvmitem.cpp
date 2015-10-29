#include "stdafx.h"
#include "mrvm.h"

//======================================================================

MRVMItem::MRVMItem(const QString& name)
{
    this->m_name = name;
}

MRVMItem::~MRVMItem()
{

}

void MRVMItem::setName(const QString& name)
{
    this->m_name = name;
}

QString MRVMItem::name() const
{
    return m_name;
}

void MRVMItem::addValue(const QString &value)
{
  this->m_values.append(value);
}

void MRVMItem::clearValues()
{
    this->m_values.clear();
}

const QList<QString>& MRVMItem::values() const
{
    return this->m_values;
}

void MRVMItem::setProperties(const QMap<QString, QString>& properties)
{
    this->m_properties = properties;
}

const QMap<QString, QString>& MRVMItem::properties() const
{
    return m_properties;
}

