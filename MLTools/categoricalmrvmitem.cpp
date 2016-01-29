#include "stdafx.h"
#include "mrvm.h"


CategoricalMRVMItem::CategoricalMRVMItem(MRVMItem::IOType iotype, const QString& name)
  :MRVMItem(iotype, name)
{

}

CategoricalMRVMItem::~CategoricalMRVMItem()
{

}

float* CategoricalMRVMItem::trainingValues(int row)
{
  float* values = new float[m_categories.count()];

  std::fill_n(values, m_categories.count() , 0.0f);
  
  int category = m_trainingValuesAsString[row].toInt();

  values[m_indexByCategory[category]] =  1.0f;
  
  return values;
}

void CategoricalMRVMItem::setTrainingValues(int row, float*& values) 
{
  ASSERT(row < m_trainingValuesAsString.count(),"Row must be less or equal to number of rows");

  int mindex = -1;
  float maxD = std::numeric_limits<float>::min();

  for(int i = 0 ; i < m_categories.size() ; i++)
    {
      float t = values[i];

      if(t > maxD)
        {
          mindex = i;
          maxD = t;
        }
    }

  if (mindex != -1)
    {
      m_trainingValuesAsString[row] = std::to_string( m_categoriesByIndex[mindex]).c_str();
    }

  delete[] values;
  values = NULL;
}

float* CategoricalMRVMItem::forecastValues(int row) 
{
  float* values = new float[m_categories.count()];

  std::fill_n(values, m_categories.count() , 0.0f);
  
  int category = m_forecastValuesAsString[row].toInt();

  values[m_indexByCategory[category]] =  1.0f;
  
  return values;
}


void CategoricalMRVMItem::setForecastValues(int row, float*& values)
{
  ASSERT(row < m_forecastValuesAsString.count(),"Row must be less or equal to number of rows");
  
  int mindex = -1;
  float maxD = std::numeric_limits<float>::min();
  
  for(int i = 0 ; i < m_categories.size() ; i++)
    {
      float t = values[i];
      
      if(t > maxD)
        {
          mindex = i;
          maxD = t;
        }
    }
  
  if (mindex != -1)
    {
      m_forecastValuesAsString[row] =  std::to_string(m_categoriesByIndex[mindex]).c_str();
    }
  
  delete[] values;
  values = NULL;
}

float* CategoricalMRVMItem::forecastUncertaintyValues(int row)
{
  ASSERT(row < m_forecastUncertaintyValuesAsString.count(),"Row must be less or equal to number of rows");

  float* values = new float[m_categories.count()];

  std::fill_n(values, m_categories.count() , std::numeric_limits<float>::max());

  int category = m_forecastValuesAsString[row].toInt();

  values[m_indexByCategory[category]] =  m_forecastUncertaintyValuesAsString[row].toFloat();

  return values;
}

void CategoricalMRVMItem::setForecastUncertaintyValues(int row, float*& values)
{
  ASSERT(row < m_forecastUncertaintyValuesAsString.count(),"Row must be less or equal to number of rows");
  int category = m_forecastValuesAsString[row].toInt();
  int index = m_indexByCategory[category];

  m_forecastUncertaintyValuesAsString[row] = std::to_string( values[index]).c_str();

  delete[] values;
  values  = NULL;
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
      else if(!xmlReader.name().compare("Categories" , Qt::CaseInsensitive))
        {
          m_categories.clear();
          m_categoriesByIndex.clear();
          m_indexByCategory.clear();
          
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
                      m_categories[key] = value.toInt();
                      m_categoriesByIndex[cc] = value.toInt();
                      m_indexByCategory[value.toInt()] = cc;
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
}

void CategoricalMRVMItem::writeXML(QXmlStreamWriter& xmlWriter)
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


  if(m_categories.size())
    {

      xmlWriter.writeStartElement("Categories");

      for(QMap<QString,int>::iterator it = m_categories.begin() ;
          it != m_categories.end() ; it++)
        {
          xmlWriter.writeStartElement("Category");

          xmlWriter.writeAttribute("Name", it.key());

          xmlWriter.writeCharacters( std::to_string( it.value()).c_str());

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

int CategoricalMRVMItem::columnCount()
{
  return m_categories.size();
}


MRVMItem::MRVMValueType CategoricalMRVMItem::valueType() const
{
  return MRVMItem::Categorical;
}

QString CategoricalMRVMItem::type() const
{
  return "CategoricalMRVMItem";
}

QMap<QString,int> CategoricalMRVMItem::categories() const
{
  return m_categories;
}


