#include <include/stdafx.h>
#include <include/mrvm.h>
#include <QTextStream>


RealMRVMItem::RealMRVMItem(MRVMItem::IOType iotype, const QString& name)
    :MRVMItem(iotype, name)
{

}

RealMRVMItem::~RealMRVMItem()
{

}

af::array RealMRVMItem::trainingValues(int valueIndex, int startRow, int length)
{
    af::array value(1,1);

    if(!m_properties["ReadFromFile"].toBool())
    {
        value = m_trainingValuesAsString[valueIndex].toDouble();
    }
    else
    {
        value = m_trainingValues[valueIndex];
    }

    return value;
}

void RealMRVMItem::setTrainingValuesAsString(const QList<QString>& trainingValues)
{
    this->m_trainingValuesAsString.clear();

    if(!m_properties["ReadFromFile"].toBool())
    {
        bool convert = false;

        for(QList<QString>::iterator it; it != trainingValues.end() ; it++)
        {
            QString valueAsString = *it;
            valueAsString.toDouble(&convert);
            ASSERT(convert , "Values in list is not a valid float");
        }
    }
    else
    {
        readWriteTrainingValuesFiles();
    }
}

af::array RealMRVMItem::forecastValues(int valueIndex, int startRow, int length)
{
    af::array value(1,1);

    if(!m_properties["ReadFromFile"].toBool())
    {
        value = m_forecastValuesAsString[valueIndex].toDouble();
    }
    else
    {
        value = m_forecastValues[valueIndex];
    }

    return value;
}

void RealMRVMItem::setForecastValues(int row, const af::array& values, const af::array& uncertainty)
{
    float* val = values.host<float>();
    float* uncert = uncertainty.host<float>();



    if(!m_properties["ReadFromFile"].toBool())
    {
        if(row >= m_forecastValuesAsString.length())
            expandListTo(m_forecastValuesAsString, row);
        if(row >= m_forecastUncertaintyValuesAsString.length())
            expandListTo(m_forecastUncertaintyValuesAsString, row);

        this->m_forecastValuesAsString[row] = QString::number(val[0]);
        this->m_forecastUncertaintyValuesAsString[row] = QString::number(uncert[0]);
    }
    else
    {
        if(row >= m_forecastValues.length())
            expandListTo(m_forecastValues, row);
        if(row >= m_forecastUncertaintyValues.length())
            expandListTo(m_forecastUncertaintyValues, row);

        this->m_forecastValues[row] = val[0];
        this->m_forecastUncertaintyValues[row] = uncert[0];
    }


    delete[] val;
    delete[] uncert;
}

void RealMRVMItem::readXML(QXmlStreamReader &xmlReader)
{
    MRVMItem::readXML(xmlReader);

    if(m_properties["ReadFromFile"].toBool())
    {
        readWriteTrainingValuesFiles();
        readWriteForecastValuesFiles();
        readWriteForecastUncertaintyValuesFiles();
    }
}

void RealMRVMItem::writeXML(QXmlStreamWriter &xmlWriter)
{
    if(m_properties["ReadFromFile"].toBool())
    {
        readWriteTrainingValuesFiles(false);
        readWriteForecastValuesFiles(false);
        readWriteForecastUncertaintyValuesFiles(false);
    }

    MRVMItem::writeXML(xmlWriter);
}

int RealMRVMItem::columnCount()
{
    return 1;
}

int RealMRVMItem::numTrainingValues() const
{
    if(m_properties["ReadFromFile"].toBool())
    {
        return m_trainingValues.count();
    }
    return m_trainingValuesAsString.count();
}

int RealMRVMItem::numForecastValues() const
{
    if(m_properties["ReadFromFile"].toBool())
    {
        return m_forecastValues.count();
    }

    return m_forecastValuesAsString.count();
}

MRVMItem::MRVMValueType RealMRVMItem::valueType() const
{
    return MRVMItem::Real;
}

QString RealMRVMItem::type() const
{
    return "RealMRVMItem";
}

void RealMRVMItem::readWriteTrainingValuesFiles(bool read)
{


    QFile file (m_trainingValuesAsString[0]);

    if(read)
    {
        if(file.open(QIODevice::ReadOnly | QIODevice::Text))
        {
            m_trainingValues.clear();
            bool ok;
            while (!file.atEnd())
            {
                QString line = file.readLine();
                float value = line.trimmed().toFloat(&ok);

                if(ok)
                {
                    m_trainingValues.append(value);
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
                tStream << m_trainingValues[i] << endl;
            }
        }
    }

    file.close();
}

void RealMRVMItem::readWriteForecastValuesFiles(bool read)
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
                while (!file.atEnd())
                {
                    QString line = file.readLine();
                    float value = line.trimmed().toFloat(&ok);

                    if(ok)
                    {
                        m_forecastValues.append(value);
                    }
                }
            }
        }
        else
        {
            if(file.open(QIODevice::WriteOnly | QIODevice::Text))
            {
                QTextStream tStream(&file);

                tStream << "Forecast Values" << endl;

                for(int i = 0 ; i < m_forecastValues.length() ; i++)
                {
                    tStream << m_forecastValues[i] << endl;
                }
            }
        }

        file.close();
    }
}

void RealMRVMItem::readWriteForecastUncertaintyValuesFiles(bool read)
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
                while (!file.atEnd())
                {
                    QString line = file.readLine();
                    float value = line.trimmed().toFloat(&ok);

                    if(ok)
                    {
                        m_forecastUncertaintyValues.append(value);
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
                    tStream << m_forecastUncertaintyValues[i] << endl;
                }
            }
        }

        file.close();
    }
}

void RealMRVMItem::expandListTo(QList<float>& list, int index)
{
    while (list.length() <= index )
    {
        list.append(0);
    }
}

void RealMRVMItem::expandListTo(QList<QString>& list, int index)
{
    while (list.length() <= index )
    {
        list.append("");
    }
}
