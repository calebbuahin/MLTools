#include "stdafx.h"
#include "mrvm.h"

RealMRVMItem::RealMRVMItem(MRVMItem::IOType iotype, const QString& name)
  :MRVMItem(iotype, name)
{

}

RealMRVMItem::~RealMRVMItem()
{

}


float* RealMRVMItem::trainingValues(int row)
{
  float* value = new float[1];
  value[0] = m_trainingValuesAsString[row].toDouble();
  return value;
}

void RealMRVMItem::setTrainingValuesAsString(const QList<QString>& trainingValues)
{
  this->m_trainingValuesAsString.clear();
  bool convert = false;

  for(QList<QString>::iterator it; it != trainingValues.end() ; it++)
    {
      QString valueAsString = *it;
      valueAsString.toDouble(&convert);
      ASSERT(convert , "Values in list is not a valid float");
    }
}


void RealMRVMItem::setTrainingValues(int row, float*& values)
{
  this->m_trainingValuesAsString[row] = QString(std::to_string(values[0]).c_str());
  delete values;
  values = NULL;
}


float* RealMRVMItem::forecastValues(int row)
{
  float* value = new float[1];
  value[0] = m_forecastValuesAsString[row].toDouble();
  return value;
}


void RealMRVMItem::setForecastValues(int index, float*& values)
{
  this->m_forecastValuesAsString[index] = QString::number(values[0]);
  delete values;
  values = NULL;
}

float* RealMRVMItem::forecastUncertaintyValues(int row)
{
  float* value = new float[1];
  value[0] = m_forecastUncertaintyValuesAsString[row].toDouble();
  return value;
}

void RealMRVMItem::setForecastUncertaintyValues(int row, float*& values)
{
  m_forecastUncertaintyValuesAsString[row] = QString::number(values[0]);
  delete values;
  values = NULL;
}

int RealMRVMItem::columnCount()
{
  return 1;
}

MRVMItem::MRVMValueType RealMRVMItem::valueType() const
{
  return MRVMItem::Real;
}

QString RealMRVMItem::type() const
{
  return "RealMRVMItem";
}
