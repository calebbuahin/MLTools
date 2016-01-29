#include "stdafx.h"
#include "mrvm.h"

//======================================================================

MRVMItem::MRVMItem(MRVMItem::IOType iotype , const QString& name)
  :m_name(name),m_iotype(iotype)
{
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
          m_properties.clear();

          while (!(xmlReader.isEndElement() && !xmlReader.name().compare("Properties" , Qt::CaseInsensitive)) && !xmlReader.hasError())
            {
              if(!xmlReader.name().compare("KeyValue", Qt::CaseInsensitive))
                {
                  QXmlStreamAttributes propAttributes =  xmlReader.attributes();
                  QString key("");
                  QString value = xmlReader.readElementText();

                  QXmlStreamAttribute attribute = *propAttributes.begin();

                  if(!attribute.name().compare("Key",Qt::CaseInsensitive))
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
}

void MRVMItem::writeXML(QXmlStreamWriter &xmlWriter)
{

  xmlWriter.writeStartElement("MRVMItem");
  xmlWriter.writeAttribute("name",m_name);
  xmlWriter.writeAttribute("type" , type());


  if(m_properties.size())
    {

      xmlWriter.writeStartElement("Properties");

      for(QMap<QString,QString>::iterator it = m_properties.begin() ;
          it != m_properties.end() ; it++)
      {
         xmlWriter.writeStartElement("KeyValue");

         xmlWriter.writeAttribute("Key", it.key());

         xmlWriter.writeCharacters(it.value());

         xmlWriter.writeEndElement();
      }

      xmlWriter.writeEndElement();
    }


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

int MRVMItem::numTrainingValues() const
{
  return m_trainingValuesAsString.count();
}

int MRVMItem::numForecastValues() const
{
  return m_forecastValuesAsString.count();
}

const QMap<QString, QString>& MRVMItem::properties() const
{
  return m_properties;
}

void MRVMItem::setProperties(const QMap<QString, QString>& properties)
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

  if(m_properties.count())
    {
      output += "  <Properties>\n";

      for(QMap<QString,QString>::const_iterator it = m_properties.begin(); it != m_properties.end() ; it++)
        {
          output += "    <KeyValue Key=\"" + it.key() + "\">" + it.value() + "/KeyValue>\n";
        }

      output += "  </Properties>\n";
    }

  if(m_trainingValuesAsString.count())
    {
      output += "  <TrainingValues>\n";

      for(QList<QString>::const_iterator it = m_trainingValuesAsString.begin(); it != m_trainingValuesAsString.end() ; it++)
        {
          output += "    <Value>" + (*it) + "/Value>\n";
        }

      output += "  </TrainingValues>\n";
    }

  if(m_forecastValuesAsString.count())
    {
      output += "  <ForecastValues>\n";

      for(QList<QString>::const_iterator it = m_forecastValuesAsString.begin(); it != m_forecastValuesAsString.end() ; it++)
        {
          output += "    <Value>" + (*it) + "/Value>\n";
        }

      output += "  </ForecastValues>\n";
    }

  if(m_forecastUncertaintyValuesAsString.count())
    {
      output += "  <ForecastUncertaintyValues>\n";

      for(QList<QString>::const_iterator it = m_forecastUncertaintyValuesAsString.begin(); it != m_forecastUncertaintyValuesAsString.end() ; it++)
        {
          output += "    <Value>" + (*it) + "/Value>\n";
        }

      output += "  </ForecastUncertaintyValues>\n";
    }

  output += "</MRVM>\n";

  return output;
}
