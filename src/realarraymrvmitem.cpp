#include <include/stdafx.h>
#include <include/mrvm.h>
#include <QTextStream>

QRegularExpression RealArrayMRVMItem::s_regex("\\s*(=>|,|\\s)\\s*");


RealArrayMRVMItem::RealArrayMRVMItem(MRVMItem::IOType iotype,const QString& name)
    :MRVMItem(iotype, name),  m_columnCount(0)
{

}

RealArrayMRVMItem::~RealArrayMRVMItem()
{

}

af::array RealArrayMRVMItem::trainingValues(int row)
{

    af::array values(1, m_columnCount);

    if(!m_properties["ReadFromFile"].toBool())
    {
        ASSERT(row < m_trainingValuesAsString.count(),"Row must be less or equal to number of rows");


        QString line = m_trainingValuesAsString[row];
        QStringList columnValues = line.split(",");

        for(int i = 0 ; i < m_columnCount; i++)
        {
            values(0,i) = columnValues[i].toFloat();
        }
    }
    else
    {
        QList<float> vals = m_trainingValues[row];

        for(int i = 0 ; i < m_columnCount; i++)
        {
            values(0,i) = vals[i];
        }

    }

    return values;
}

void RealArrayMRVMItem::setTrainingValuesAsString(const QList<QString> &trainingValues)
{
    MRVMItem::setTrainingValuesAsString(trainingValues);

    if(!m_properties["ReadFromFile"].toBool())
    {
        QStringList line = trainingValues[0].split(",");
        m_columnCount = line.size();
    }
    else
    {
        readWriteTrainingValuesFiles();
    }
}

af::array RealArrayMRVMItem::forecastValues(int row)
{

    af::array values(1, m_columnCount);

    if(!m_properties["ReadFromFile"].toBool())
    {
        ASSERT(row < m_forecastValuesAsString.count(),"Row must be less or equal to number of rows");

        QString line = m_forecastValuesAsString[row];
        QStringList columnValues = line.split(",");

        for(int i = 0 ; i < m_columnCount; i++)
        {
            values(0,i) = columnValues[i].toFloat();
        }
    }
    else
    {
        QList<float> vals = m_trainingValues[row];

        for(int i = 0 ; i < m_columnCount; i++)
        {
            values(0,i) = vals[i];
        }
    }

    return values;
}

void RealArrayMRVMItem::setForecastValues(int row, const af::array& values, const af::array& uncertainty)
{
    float* val = values.host<float>();
    float* uncert = uncertainty.host<float>();

    if(!m_properties["ReadFromFile"].toBool())
    {
        if(row >= m_forecastValuesAsString.length())
            expandListTo(m_forecastValuesAsString, row);
        if(row >= m_forecastUncertaintyValuesAsString.length())
            expandListTo(m_forecastUncertaintyValuesAsString, row);

        QString line =  QString::number(val[0]);
        QString lineUncert =  QString::number(uncert[0]);

        if(m_columnCount > 1)
        {
            for(int i = 1 ; i < m_columnCount; i++)
            {
                line = line + "," + QString::number(val[i]);
                lineUncert = lineUncert + "," + QString::number(uncert[i]);
            }
        }

        m_forecastValuesAsString[row] = line;
        m_forecastUncertaintyValuesAsString[row] = lineUncert;
    }
    else
    {
        if(row >= m_forecastValues.length())
            expandListTo(m_forecastValues,row);

        QList<float> rowV = m_forecastValues[row];

        for(int i = 0 ; i < m_columnCount; i++)
        {
            rowV[i] = val[i];
        }

        if(row >= m_forecastUncertaintyValues.length())
            expandListTo(m_forecastUncertaintyValues,row);

        QList<float> rowU = m_forecastUncertaintyValues[row];

        for(int i = 0 ; i < m_columnCount; i++)
        {
            rowU[i] = uncert[i];
        }
    }

    delete[] val;
    delete[] uncert;
}

void RealArrayMRVMItem::readXML(QXmlStreamReader &xmlReader)
{

    MRVMItem::readXML(xmlReader);

    if(m_properties["ReadFromFile"].toBool())
    {
        readWriteTrainingValuesFiles();
        readWriteForecastValuesFiles();
        readWriteForecastUncertaintyValuesFiles();
    }
    else
    {
        QStringList line = m_trainingValuesAsString[0].split(",");
        m_columnCount = line.size();
    }
}

void RealArrayMRVMItem::writeXML(QXmlStreamWriter &xmlWriter)
{
    if(m_properties["ReadFromFile"].toBool())
    {
        readWriteTrainingValuesFiles(false);
        readWriteForecastValuesFiles(false);
        readWriteForecastUncertaintyValuesFiles(false);
    }

    MRVMItem::writeXML(xmlWriter);
}

int RealArrayMRVMItem::columnCount()
{
    return m_columnCount;
}

int RealArrayMRVMItem::numTrainingValues() const
{
    if(m_properties["ReadFromFile"].toBool())
    {
        return m_trainingValues.count();
    }
    return m_trainingValuesAsString.count();
}

int RealArrayMRVMItem::numForecastValues() const
{
    if(m_properties["ReadFromFile"].toBool())
    {
        return m_forecastValues.count();
    }

    return m_forecastValuesAsString.count();
}

MRVMItem::MRVMValueType RealArrayMRVMItem::valueType() const
{
    return MRVMItem::Real;
}

QString RealArrayMRVMItem::type() const
{
    return "RealArrayMRVMItem";
}

void RealArrayMRVMItem::readWriteTrainingValuesFiles(bool read)
{
    QFile file (m_trainingValuesAsString[0]);

    if(read)
    {
        if(file.open(QIODevice::ReadOnly | QIODevice::Text))
        {
            m_trainingValues.clear();
            bool ok;
            m_columnCount = 0;

            while (!file.atEnd())
            {
                QString line = file.readLine();
                QStringList  lineList = line.split(s_regex);

                QList<float> mvals;

                for(int i = 0 ; i < lineList.count() ; i++)
                {
                    float value = lineList[i].toFloat(&ok);

                    if(ok)
                    {
                        mvals.append(value);
                    }
                    else
                    {
                        break;
                    }
                }

                if(mvals.count())
                {
                    m_trainingValues.append(mvals);

                    if(!m_columnCount)
                    {
                        m_columnCount = mvals.count();
                    }
                    else
                    {
                        ASSERT(m_columnCount == mvals.count() , "Column count mismatches");
                    }
                }
            }
        }
    }
    else
    {
        if(file.open(QIODevice::WriteOnly | QIODevice::Text))
        {
            QTextStream tStream(&file);
            tStream << "Training Values" << endl;

            for(int i = 0 ; i < m_trainingValues.length() ; i++)
            {
                QList<float> row = m_trainingValues[i];

                for(int c = 0 ; c < row.length() ; c++)
                {
                    if(c)
                    {
                        tStream << "," << m_trainingValues[i][c];
                    }
                    else
                    {
                        tStream << m_trainingValues[i][c];
                    }

                }

                tStream << endl;
            }
        }
    }

    file.close();
}

void RealArrayMRVMItem::readWriteForecastValuesFiles(bool read)
{

    if(m_forecastValuesAsString.length())
    {
        QFile file (m_forecastValuesAsString[0]);

        if(read)
        {
            if(file.open(QIODevice::ReadOnly | QIODevice::Text))
            {
                m_forecastValues.clear();
                bool ok;
                m_columnCount = 0;

                while (!file.atEnd())
                {
                    QString line = file.readLine();
                    QStringList  lineList = line.split(s_regex);

                    QList<float> mvals;

                    for(int i = 0 ; i < lineList.count() ; i++)
                    {
                        float value = lineList[i].toFloat(&ok);

                        if(ok)
                        {
                            mvals.append(value);
                        }
                    }

                    if(mvals.count())
                    {
                        m_forecastValues.append(mvals);

                        if(!m_columnCount)
                        {
                            m_columnCount = mvals.count();
                        }
                        else
                        {
                            ASSERT(m_columnCount == mvals.count() , "Column count mismatches");
                        }
                    }
                }
            }
        }
        else
        {
            if(file.open(QIODevice::WriteOnly | QIODevice::Text))
            {
                QTextStream tStream(&file);
                tStream << "Training Values" << endl;

                for(int i = 0 ; i < m_forecastValues.length() ; i++)
                {
                    QList<float> row = m_forecastValues[i];

                    for(int c = 0 ; c < row.length() ; c++)
                    {
                        if(c)
                        {
                            tStream << "," << m_forecastValues[i][c];
                        }
                        else
                        {
                            tStream << m_forecastValues[i][c];
                        }

                    }

                    tStream << endl;
                }
            }
        }


        file.close();
    }
}

void RealArrayMRVMItem::readWriteForecastUncertaintyValuesFiles(bool read)
{

    if(m_forecastUncertaintyValuesAsString.length())
    {
        QFile file (m_forecastUncertaintyValuesAsString[0]);

        if(read)
        {
            if(file.open(QIODevice::ReadOnly | QIODevice::Text))
            {
                m_forecastUncertaintyValues.clear();
                bool ok;
                m_columnCount = 0;

                while (!file.atEnd())
                {
                    QString line = file.readLine();
                    QStringList  lineList = line.split(s_regex);

                    QList<float> mvals;

                    for(int i = 0 ; i < lineList.count() ; i++)
                    {
                        float value = lineList[i].toFloat(&ok);

                        if(ok)
                        {
                            mvals.append(value);
                        }
                    }

                    if(mvals.count())
                    {
                        m_forecastUncertaintyValues.append(mvals);

                        if(!m_columnCount)
                        {
                            m_columnCount = mvals.count();
                        }
                        else
                        {
                            ASSERT(m_columnCount == mvals.count() , "Column count mismatches");
                        }
                    }
                }
            }
        }
        else
        {
            if(file.open(QIODevice::WriteOnly | QIODevice::Text))
            {
                QTextStream tStream(&file);
                tStream << "Forecast Uncertainty Values" << endl;

                for(int i = 0 ; i < m_forecastUncertaintyValues.length() ; i++)
                {
                    QList<float> row = m_forecastUncertaintyValues[i];

                    for(int c = 0 ; c < row.length() ; c++)
                    {
                        if(c)
                        {
                            tStream << "," << m_forecastUncertaintyValues[i][c];
                        }
                        else
                        {
                            tStream << m_forecastUncertaintyValues[i][c];
                        }

                    }

                    tStream << endl;
                }
            }
        }

        file.close();
    }
}

void RealArrayMRVMItem::expandListTo(QList<QList<float> >& list, int index)
{
    while (list.length() < index )
    {
        QList<float> values;

        for(int i = 0 ; i < m_columnCount ; i++)
            values.append(0);

        list.append(values);
    }
}

void RealArrayMRVMItem::expandListTo(QList<QString>& list, int index)
{
    while (list.length() <= index )
    {
        list.append("");
    }
}
