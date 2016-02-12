#include <headers/stdafx.h>
#include <headers/mrvm.h>

//======================================================================

MRVMItem::MRVMItem(MRVMItem::IOType iotype , const QString& name)
  :m_iotype(iotype), m_name(name)
{
  m_properties["ReadFromFile"] = false;
}

MRVMItem::~MRVMItem()
{

}

QString MRVMItem::name() const
{
  return m_name;
}

void MRVMItem::setName(const QString& name)
{
  this->m_name = name;
}

void MRVMItem::clearAllValues()
{
  this->m_trainingValuesAsString.clear();
  this->m_forecastValuesAsString.clear();
  this->m_forecastUncertaintyValuesAsString.clear();
}

const QList<QString>& MRVMItem::trainingValuesAsString() const
{
  return this->m_trainingValuesAsString;
}

void MRVMItem::setTrainingValuesAsString(const QList<QString>& trainingValues)
{
  this->m_trainingValuesAsString = trainingValues;
  m_numTrainingValues = trainingValues.length();
}

const QList<QString>& MRVMItem::forecastValuesAsString() const
{
  return this->m_forecastValuesAsString;
}

void MRVMItem::setForecastValuesAsString(const QList<QString>& forecastValues)
{
  this->m_forecastValuesAsString = forecastValues;
  this->m_forecastUncertaintyValuesAsString.clear();

  for(int i = 0 ; i < this->m_forecastValuesAsString.count() ; i++)
    {
      this->m_forecastUncertaintyValuesAsString.append(QString());
    }

  m_numForecastValues = forecastValues.length() ;
}

const QList<QString>& MRVMItem::forecastUncertaintyValuesAsString() const
{
  return this->m_forecastUncertaintyValuesAsString;
}

void MRVMItem::setForecastUncertaintyValueAsString(const QList<QString>& forecastUncertaintyValuesAsString)
{
  this->m_forecastUncertaintyValuesAsString = forecastUncertaintyValuesAsString;
}

void MRVMItem::readXML(QXmlStreamReader & xmlReader)
{
  while (!(xmlReader.isEndElement() && !xmlReader.name().compare("MRVMItem", Qt::CaseInsensitive)) && !xmlReader.hasError())
    {
      if(!xmlReader.name().compare("Properties" , Qt::CaseInsensitive))
        {
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

   m_numTrainingValues = m_trainingValuesAsString.length();
   m_numForecastValues = m_forecastUncertaintyValuesAsString.length();
}

void MRVMItem::writeXML(QXmlStreamWriter &xmlWriter)
{

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


  xmlWriter.writeStartElement("TrainingValues");

  if(m_trainingValuesAsString.count())
    {

      for(QList<QString>::iterator it = m_trainingValuesAsString.begin() ; it != m_trainingValuesAsString.end() ; it++)
        {
          xmlWriter.writeTextElement("Value",*it);
        }
    }

  xmlWriter.writeEndElement();


  xmlWriter.writeStartElement("ForecastValues");

  if(m_forecastValuesAsString.count())
    {
      for(QList<QString>::iterator it = m_forecastValuesAsString.begin(); it != m_forecastValuesAsString.end() ; it++)
        {
          xmlWriter.writeTextElement("Value",*it);
        }
    }

  xmlWriter.writeEndElement();


  xmlWriter.writeStartElement("ForecastUncertaintyValues");

  if(m_forecastUncertaintyValuesAsString.count())
    {

      for(QList<QString>::iterator it = m_forecastUncertaintyValuesAsString.begin() ; it != m_forecastUncertaintyValuesAsString.end() ; it++)
        {
          xmlWriter.writeTextElement("Value",*it);
        }

    }

  xmlWriter.writeEndElement();

  xmlWriter.writeEndElement();
}

int MRVMItem::numTrainingValues() const
{
  return m_numTrainingValues;
}

int MRVMItem::numForecastValues() const
{
  return m_numForecastValues;
}

int MRVMItem::numRowsPerTrainingValue() const
{
  return 1;
}

int MRVMItem::numRowsPerForecastValue() const
{
  return 1;
}

const QMap<QString, QVariant>& MRVMItem::properties() const
{
  return m_properties;
}

void MRVMItem::setProperties(const QMap<QString, QVariant>& properties)
{
  this->m_properties = properties;
}

MRVMItem::IOType MRVMItem::ioType() const
{
  return m_iotype;
}

QString MRVMItem::toString() const
{

  QString output;

  output += "<MRVM name=\"" + m_name + "\" type=\"" + type() + "\">\n";

  output += "  <Properties>\n";

  if(m_properties.count())
    {
      for(QMap<QString,QVariant>::const_iterator it = m_properties.begin(); it != m_properties.end() ; it++)
        {
          output += "    <Property Key=\"" + it.key() + "\">" + it.value().toString() + "</Property>\n";
        }
    }

  output += "  </Properties>\n";



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
