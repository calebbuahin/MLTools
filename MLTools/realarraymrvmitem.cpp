#include "stdafx.h"
#include "mrvm.h"

RealArrayMRVMItem::RealArrayMRVMItem(MRVMItem::IOType iotype,const QString& name)
  :MRVMItem(iotype, name),  m_columnCount(0)
{

}

RealArrayMRVMItem::~RealArrayMRVMItem()
{

}


float* RealArrayMRVMItem::trainingValues(int row)
{
  ASSERT(row < m_trainingValuesAsString.count(),"Row must be less or equal to number of rows");

  QString line = m_trainingValuesAsString[row];
  QStringList columnValues = line.split(",");

  float* values = new float[m_columnCount];

  for(int i = 0 ; i < m_columnCount; i++)
    {
      values[i] = columnValues[i].toFloat();
    }

  return values;
}

void RealArrayMRVMItem::setTrainingValuesAsString(const QList<QString> &trainingValues)
{
  QStringList line = trainingValues[0].split(",");

  m_columnCount = line.size();

  MRVMItem::setTrainingValuesAsString(trainingValues);
}

void RealArrayMRVMItem::setTrainingValues(int row, float *&values)
{
  ASSERT(row < m_trainingValuesAsString.count(),"Row must be less or equal to number of rows");

  QString line =  QString::number(values[0]);

  if(m_columnCount > 1)
    {
      for(int i = 1 ; i < m_columnCount; i++)
        {
          line = line + "," + QString::number(values[i]);
        }
    }

  m_trainingValuesAsString[row] = line;

  delete[] values;
  values = NULL;
}


float* RealArrayMRVMItem::forecastValues(int row)
{
  ASSERT(row < m_forecastValuesAsString.count(),"Row must be less or equal to number of rows");

  QString line = m_forecastValuesAsString[row];
  QStringList columnValues = line.split(",");

  float* values = new float[m_columnCount];

  for(int i = 0 ; i < m_columnCount; i++)
    {
      values[i] = columnValues[i].toFloat();
    }

  return values;
}

void RealArrayMRVMItem::setForecastValues(int row, float *&values)
{
  ASSERT(row < m_forecastValuesAsString.count(),"Row must be less or equal to number of rows");

  QString line =  QString::number(values[0]);

  if(m_columnCount > 1)
    {
      for(int i = 1 ; i < m_columnCount; i++)
        {
          line = line + "," + QString::number(values[i]);
        }
    }

  m_forecastValuesAsString[row] = line;

  delete[] values;
  values = NULL;
}

float* RealArrayMRVMItem::forecastUncertaintyValues(int row)
{
  ASSERT(row < m_forecastUncertaintyValuesAsString.count(),"Row must be less or equal to number of rows");

  QString line = m_forecastUncertaintyValuesAsString[row];
  QStringList columnValues = line.split(",");

  float* values = new float[m_columnCount];

  for(int i = 0 ; i < m_columnCount; i++)
    {
      values[i] = columnValues[i].toFloat();
    }

  return values;
}

void RealArrayMRVMItem::setForecastUncertaintyValues(int row, float*& values)
{
  ASSERT(row < m_forecastUncertaintyValuesAsString.count(),"Row must be less or equal to number of rows");

  QString line =  QString::number(values[0]);

  if(m_columnCount > 1)
    {
      for(int i = 1 ; i < m_columnCount; i++)
        {
          line = line + "," + QString::number(values[i]);
        }
    }

  m_forecastUncertaintyValuesAsString[row] = line;

  delete[] values;
  values = NULL;
}

void RealArrayMRVMItem::readXML(QXmlStreamReader &xmlReader)
{
  MRVMItem::readXML(xmlReader);
  QStringList line = m_trainingValuesAsString[0].split(",");
  m_columnCount = line.size();
}

int RealArrayMRVMItem::columnCount()
{
  return m_columnCount;
}

MRVMItem::MRVMValueType RealArrayMRVMItem::valueType() const
{
  return MRVMItem::Real;
}

QString RealArrayMRVMItem::type() const
{
  return "RealArrayMRVMItem";
}
