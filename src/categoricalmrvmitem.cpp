#include <include/stdafx.h>
#include <include/mrvm.h>
#include <QTextStream>

CategoricalMRVMItem::CategoricalMRVMItem(MRVMItem::IOType iotype, const QString& name)
    :MRVMItem(iotype, name), maxCValue(100.0f),minCValue(0.0000001f)
{

}

CategoricalMRVMItem::~CategoricalMRVMItem()
{

}

QMap<QString,int> CategoricalMRVMItem::getCategories() const
{
    return m_classbycategory;
}

void CategoricalMRVMItem::setCategories(const QMap<QString, int> &categories)
{
    this->m_classbycategory = categories;
}

af::array CategoricalMRVMItem::trainingValues(int valueIndex, int startRow, int length)
{
    af::array values(1, m_classbycategory.count());

    values =  minCValue;

    if(!m_properties["ReadFromFile"].toBool())
    {
        int category = m_trainingValuesAsString[valueIndex].toInt();
        values(0,m_indexbyclass[category]) =  maxCValue;
    }
    else
    {
        int category = m_trainingValues[valueIndex];
        values(0,m_indexbyclass[category]) =  maxCValue;
    }

    return values;
}

af::array CategoricalMRVMItem::forecastValues(int valueIndex, int startRow, int length)
{
    af::array values(1,m_classbycategory.count());

    values =  minCValue;

    if(!m_properties["ReadFromFile"].toBool())
    {
        int category = m_forecastValuesAsString[valueIndex].toInt();
        values(0, m_indexbyclass[category]) =  maxCValue;
    }
    else
    {
        int category = m_forecastValues[valueIndex];
        values(0, m_indexbyclass[category]) =  maxCValue;
    }

    return values;
}

void CategoricalMRVMItem::setForecastValues(int row, const af::array& values, const af::array& uncertainty)
{
    int mindex = -1;
    float maxD = std::numeric_limits<float>::min();

    float* val = values.host<float>();
    float* uncert = uncertainty.host<float>();

    for(int i = 0 ; i < m_indexbyclass.size() ; i++)
    {
        float t = val[i];

        if(t > maxD)
        {
            mindex = i;
            maxD = t;
        }
    }

    if (mindex != -1)
    {
        if(!m_properties["ReadFromFile"].toBool())
        {
            if(row >= m_forecastValuesAsString.length())
                expandListTo(m_forecastValuesAsString, row);
            if(row >= m_forecastUncertaintyValuesAsString.length())
                expandListTo(m_forecastUncertaintyValuesAsString, row);

            m_forecastValuesAsString[row] =  QString::number(m_classbyindex[mindex]);
            m_forecastUncertaintyValuesAsString[row] = QString::number(uncert[mindex]);
        }
        else
        {
            if(row >= m_forecastValues.length())
                expandListTo(m_forecastValues, row);

            m_forecastValues[row] = m_classbyindex[mindex];

            if(row >= m_forecastUncertaintyValues.length())
                expandListTo(m_forecastUncertaintyValues, row);

            m_forecastUncertaintyValues[row] = uncert[mindex];
        }
    }

    delete[] val;
}

void CategoricalMRVMItem::readXML(QXmlStreamReader & xmlReader)
{
    while (!(xmlReader.isEndElement() && !xmlReader.name().compare("MRVMItem", Qt::CaseInsensitive)) && !xmlReader.hasError())
    {
        if(!xmlReader.name().compare("Properties" , Qt::CaseInsensitive))
        {
            m_properties.clear();

            while (!(xmlReader.isEndElement() && !xmlReader.name().compare("Properties" , Qt::CaseInsensitive)) && !xmlReader.hasError())
            {
                if(!xmlReader.name().compare("Property", Qt::CaseInsensitive))
                {
                    QXmlStreamAttributes propAttributes =  xmlReader.attributes();
                    QString key("");
                    QString value = xmlReader.readElementText();

                    QXmlStreamAttribute attribute = *propAttributes.begin();

                    if(!attribute.name().compare("Name",Qt::CaseInsensitive))
                    {
                        key = attribute.value().toString();
                    }

                    if(!key.isNull()  && !key.isEmpty() &&  ! value.isNull() && !value.isEmpty())
                    {
                        m_properties[key] = value;
                    }
                }
                xmlReader.readNext() ;
            }
        }
        else if(!xmlReader.name().compare("Categories" , Qt::CaseInsensitive))
        {
            m_classbycategory.clear();
            m_indexbyclass.clear();
            m_classbyindex.clear();

            int cc = 0;

            while (!(xmlReader.isEndElement() && !xmlReader.name().compare("Categories" , Qt::CaseInsensitive)) && !xmlReader.hasError())
            {
                if(!xmlReader.name().compare("Category", Qt::CaseInsensitive))
                {
                    QXmlStreamAttributes propAttributes =  xmlReader.attributes();
                    QString key("");
                    QString value = xmlReader.readElementText();

                    QXmlStreamAttribute attribute = *propAttributes.begin();

                    if(!attribute.name().compare("Name",Qt::CaseInsensitive))
                    {
                        key = attribute.value().toString();
                    }

                    if(!key.isNull()  && !key.isEmpty() &&  ! value.isNull() && !value.isEmpty())
                    {
                        m_classbycategory[key] = value.toInt();
                        m_categorybyclass[value.toInt()] = key;
                        m_classbyindex[cc] = value.toInt();
                        m_indexbyclass[value.toInt()] = cc;
                        cc++;
                    }
                }
                xmlReader.readNext() ;
            }
        }
        else if(!xmlReader.name().compare("TrainingValues" , Qt::CaseInsensitive))
        {
            this->m_trainingValuesAsString.clear();

            while (!(xmlReader.isEndElement() && !xmlReader.name().compare("TrainingValues" , Qt::CaseInsensitive)) && !xmlReader.hasError())
            {
                if(!xmlReader.name().compare("Value", Qt::CaseInsensitive))
                {
                    this->m_trainingValuesAsString.append(xmlReader.readElementText());
                }

                xmlReader.readNext();
            }

        }
        else if(!xmlReader.name().compare("ForecastValues" , Qt::CaseInsensitive))
        {
            this->m_forecastValuesAsString.clear();

            while (!(xmlReader.isEndElement() && !xmlReader.name().compare("ForecastValues" , Qt::CaseInsensitive)) && !xmlReader.hasError())
            {
                if(!xmlReader.name().compare("Value", Qt::CaseInsensitive))
                {
                    this->m_forecastValuesAsString.append(xmlReader.readElementText());
                }

                xmlReader.readNext();
            }
        }
        else if(!xmlReader.name().compare("ForecastUncertaintyValues" , Qt::CaseInsensitive))
        {
            this->m_forecastUncertaintyValuesAsString.clear();

            while (!(xmlReader.isEndElement() && !xmlReader.name().compare("ForecastUncertaintyValues" , Qt::CaseInsensitive)) && !xmlReader.hasError())
            {
                if(!xmlReader.name().compare("Value", Qt::CaseInsensitive))
                {
                    this->m_forecastUncertaintyValuesAsString.append(xmlReader.readElementText());
                }

                xmlReader.readNext();
            }
        }

        xmlReader.readNext();
    }


    if(m_properties["ReadFromFile"].toBool())
    {
        readWriteTrainingValuesFiles();
        readWriteForecastValuesFiles();
        readWriteForecastUncertaintyValuesFiles();
    }

}

void CategoricalMRVMItem::writeXML(QXmlStreamWriter& xmlWriter)
{

    if(m_properties["ReadFromFile"].toBool())
    {
        readWriteTrainingValuesFiles(false);
        readWriteForecastValuesFiles(false);
        readWriteForecastUncertaintyValuesFiles(false);
    }

    xmlWriter.writeStartElement("MRVMItem");
    xmlWriter.writeAttribute("name",m_name);
    xmlWriter.writeAttribute("type" , type());


    xmlWriter.writeStartElement("Properties");

    if(m_properties.size())
    {
        for(QMap<QString,QVariant>::iterator it = m_properties.begin() ;
            it != m_properties.end() ; it++)
        {
            xmlWriter.writeStartElement("Property");

            xmlWriter.writeAttribute("Name", it.key());

            xmlWriter.writeCharacters(it.value().toString());

            xmlWriter.writeEndElement();
        }
    }

    xmlWriter.writeEndElement();



    xmlWriter.writeStartElement("Categories");

    if(m_classbycategory.size())
    {
        for(QMap<QString,int>::iterator it = m_classbycategory.begin() ;
            it != m_classbycategory.end() ; it++)
        {
            xmlWriter.writeStartElement("Category");

            xmlWriter.writeAttribute("Name", it.key());

            xmlWriter.writeCharacters( std::to_string( it.value()).c_str());

            xmlWriter.writeEndElement();
        }
    }

    xmlWriter.writeEndElement();


    if(m_trainingValuesAsString.count())
    {
        xmlWriter.writeStartElement("TrainingValues");

        for(QList<QString>::iterator it = m_trainingValuesAsString.begin() ; it != m_trainingValuesAsString.end() ; it++)
        {
            xmlWriter.writeTextElement("Value",*it);
        }

        xmlWriter.writeEndElement();
    }

    if(m_forecastValuesAsString.count())
    {
        xmlWriter.writeStartElement("ForecastValues");

        for(QList<QString>::iterator it = m_forecastValuesAsString.begin(); it != m_forecastValuesAsString.end() ; it++)
        {
            xmlWriter.writeTextElement("Value",*it);
        }

        xmlWriter.writeEndElement();
    }

    if(m_forecastUncertaintyValuesAsString.count())
    {
        xmlWriter.writeStartElement("ForecastUncertaintyValues");

        for(QList<QString>::iterator it = m_forecastUncertaintyValuesAsString.begin() ; it != m_forecastUncertaintyValuesAsString.end() ; it++)
        {
            xmlWriter.writeTextElement("Value",*it);
        }

        xmlWriter.writeEndElement();
    }

    xmlWriter.writeEndElement();
}

int CategoricalMRVMItem::columnCount()
{
    return m_classbycategory.size();
}

int CategoricalMRVMItem::numTrainingValues() const
{
    if(m_properties["ReadFromFile"].toBool())
    {
        return m_trainingValues.count();
    }
    return m_trainingValuesAsString.count();
}

int CategoricalMRVMItem::numForecastValues() const
{
    if(m_properties["ReadFromFile"].toBool())
    {
        return m_forecastValues.count();
    }

    return m_forecastValuesAsString.count();
}

MRVMItem::MRVMValueType CategoricalMRVMItem::valueType() const
{
    return MRVMItem::Categorical;
}

QString CategoricalMRVMItem::type() const
{
    return "CategoricalMRVMItem";
}

QString CategoricalMRVMItem::toString() const
{

    QString output;

    output += "<MRVM name=\"" + m_name + "\" type=\"" + type() + "\">\n";

    output += "  <Properties>\n";

    if(m_properties.count())
    {
        for(QMap<QString,QVariant>::const_iterator it = m_properties.begin(); it != m_properties.end() ; it++)
        {
            output += "    <Property Name=\"" + it.key() + "\">" + it.value().toString() + "</Property>\n";
        }
    }

    output += "  </Properties>\n";

    output += "  <Categories>\n";

    if(m_classbycategory.count())
    {
        for(QMap<QString,int>::const_iterator it = m_classbycategory.begin(); it != m_classbycategory.end() ; it++)
        {
            output += "    <Category Key=\"" + it.key() + "\">" + QString::number(it.value()) + "</Category>\n";
        }
    }
    output += "  </Categories>\n";

    output += "  <TrainingValues>\n";

    if(m_trainingValuesAsString.count())
    {
        for(QList<QString>::const_iterator it = m_trainingValuesAsString.begin(); it != m_trainingValuesAsString.end() ; it++)
        {
            output += "    <Value>" + (*it) + "</Value>\n";
        }
    }

    output += "  </TrainingValues>\n";



    output += "  <ForecastValues>\n";

    if(m_forecastValuesAsString.count())
    {
        for(QList<QString>::const_iterator it = m_forecastValuesAsString.begin(); it != m_forecastValuesAsString.end() ; it++)
        {
            output += "    <Value>" + (*it) + "</Value>\n";
        }
    }

    output += "  </ForecastValues>\n";


    output += "  <ForecastUncertaintyValues>\n";

    if(m_forecastUncertaintyValuesAsString.count())
    {

        for(QList<QString>::const_iterator it = m_forecastUncertaintyValuesAsString.begin(); it != m_forecastUncertaintyValuesAsString.end() ; it++)
        {
            output += "    <Value>" + (*it) + "</Value>\n";
        }
    }

    output += "  </ForecastUncertaintyValues>\n";

    output += "</MRVM>\n";

    return output;

}

QMap<QString,int> CategoricalMRVMItem::categories() const
{
    return m_classbycategory;
}

void CategoricalMRVMItem::readWriteTrainingValuesFiles(bool read)
{

    QFile file (m_trainingValuesAsString[0]);

    if(read)
    {
        if(file.open(QIODevice::ReadOnly | QIODevice::Text))
        {
            m_trainingValues.clear();

            while (!file.atEnd())
            {
                QString line = file.readLine();

                if(m_classbycategory.contains(line.trimmed()))
                {
                    int classV  = m_classbycategory[line.trimmed()];
                    m_trainingValues.append(classV);
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
                tStream << m_categorybyclass[m_trainingValues[i]] << endl;
            }
        }
    }

    file.close();
}

void CategoricalMRVMItem::readWriteForecastValuesFiles(bool read)
{
    if(m_forecastValuesAsString.length())
    {
        QFile file (m_forecastValuesAsString[0]);

        if(read)
        {
            if(file.open(QIODevice::ReadOnly | QIODevice::Text))
            {
                m_forecastValues.clear();
                while (!file.atEnd())
                {
                    QString line = file.readLine();

                    if(m_classbycategory.contains(line.trimmed()))
                    {
                        int classV = m_classbycategory[line.trimmed()];
                        m_forecastValues.append(classV);
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
                    tStream << m_categorybyclass[m_forecastValues[i]] << endl;
                }
            }
        }

        file.close();
    }
}

void CategoricalMRVMItem::readWriteForecastUncertaintyValuesFiles(bool read)
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
                    float value = line.toFloat(&ok);

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

void CategoricalMRVMItem::expandListTo(QList<int>& list, int index)
{
    while (list.length() <= index )
    {
        list.append(0);
    }
}

void CategoricalMRVMItem::expandListTo(QList<float>& list, int index)
{
    while (list.length() <= index )
    {
        list.append(0);
    }
}

void CategoricalMRVMItem::expandListTo(QList<QString>& list, int index)
{
    while (list.length() <= index )
    {
        list.append("");
    }
}
