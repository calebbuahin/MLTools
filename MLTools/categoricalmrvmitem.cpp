#include "stdafx.h"
#include "mrvm.h"


CategoricalMRVMItem::CategoricalMRVMItem(const QString& name)
    :MRVMItem(name)
{

}

CategoricalMRVMItem::~CategoricalMRVMItem()
{

}

double* CategoricalMRVMItem::values(int index)
{
    if(columnCount() > 0 && index < m_values.length())
    {
        double* values = new double[columnCount()];

        int cindex = m_values[index].toInt();

        int count = 0;

        for(QMap<int,QString>::iterator  i = m_categories.begin() ; i != m_categories.end() ; i++)
        {
            if(i.key() == cindex)
            {
                values[count] = 1000.0;
            }
            else
            {
                values[count] = 0;
            }

            count++;
        }

        return values;

    }

    return NULL;
}

void CategoricalMRVMItem::setValues(int index, double*& values)
{
    if(columnCount() > 0 && index < m_values.length())
    {
        int maxKey = -9999;
        double maxValue = std::numeric_limits<double>::min();

        int count = 0;

        for(QMap<int,QString>::iterator  i = m_categories.begin() ; i != m_categories.end() ; i++)
        {
            if(values[count] > maxValue)
            {
                maxKey = i.key();
                maxValue = values[count];
            }
            count++;
        }

        m_values[index] = QString("%").arg(maxKey);
    }

    delete[] values;
    values = NULL;
}

int CategoricalMRVMItem::columnCount()
{
    return m_categories.size();
}

MRVMItem::MRVMValueType CategoricalMRVMItem::mRVMValueType() const
{
    return MRVMItem::Categorical;
}

QString CategoricalMRVMItem::type() const
{
    return "CategoricalMRVMItem";
}

const QMap<int,QString>& CategoricalMRVMItem::categories() const
{
    return m_categories;
}

void CategoricalMRVMItem::setCategories(const QMap<int,QString>& categories)
{
    this->m_categories = categories;
}

